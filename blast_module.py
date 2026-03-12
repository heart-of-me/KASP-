"""
NCBI 远程 BLAST 模块（优化版 v2）
=================================
功能：
  1. NCBI 远程 BLAST：引物/序列特异性验证（小麦 A/B/D 亚基因组分辨）
  2. 同源序列提取 + 多序列比对：找出亚基因组间差异位点（SNP/InDel）
  3. 特异性引物设计引导：基于差异位点生成引物候选区域图谱
  4. 本地 BLAST 数据库支持（可选）：用户可自建 FASTA 数据库加速

依赖（必须）：pip install biopython
可选：NCBI BLAST+ 命令行工具（本地比对模式）

v2 相对 v1 的改进：
  - 新增 BLAST 结果缓存，避免重复查询
  - 新增本地 BLAST 数据库支持（makeblastdb + blastn 命令行）
  - 新增同源序列多序列比对功能（提取 A/B/D 同源拷贝，Pairwise → 差异位点表）
  - 新增引物对设计（Forward + Reverse）在特异区段内的完整候选生成
  - 优化特异性评分：考虑 3'端匹配、错配位置、命中分散度等
  - 增加重试 + 错误恢复机制
"""

import os
import time
import re
import subprocess
import hashlib
import tempfile
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass, field

try:
    import importlib
    _bio_blast = importlib.util.find_spec("Bio")
    BIOPYTHON_AVAILABLE = _bio_blast is not None
except Exception:
    BIOPYTHON_AVAILABLE = False

def _import_bio():
    """按需导入 biopython 模块"""
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez
    return NCBIWWW, NCBIXML, Entrez


# ============================================================
# 数据类
# ============================================================

@dataclass
class HSPInfo:
    """单个 High-scoring Segment Pair"""
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    identity: float        # 0–100
    identity_count: int    # 实际匹配碱基数
    align_length: int
    mismatches: int
    gaps: int
    score: float
    e_value: float
    query_seq: str = ""
    subject_seq: str = ""
    midline: str = ""


@dataclass
class BlastHit:
    """单个 BLAST 命中"""
    title: str
    accession: str
    subject_length: int
    best_score: float
    best_e_value: float
    best_identity: float
    query_coverage: float
    subgenome: str
    chromosome: str
    hsps: List[HSPInfo] = field(default_factory=list)
    three_prime_match: bool = True
    mismatch_positions: List[int] = field(default_factory=list)


@dataclass
class BlastResult:
    """完整 BLAST 搜索结果"""
    query: str
    query_length: int
    database: str
    total_hits: int
    hits: List[BlastHit] = field(default_factory=list)
    is_specific: bool = True
    specificity_score: float = 100.0
    subgenome_hits: Dict[str, int] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    suggestions: List[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0


@dataclass
class SNPSite:
    """A/B/D 同源序列间的差异位点"""
    position: int
    ref_base: str
    alt_bases: Dict[str, str] = field(default_factory=dict)
    is_indel: bool = False
    flanking_5: str = ""
    flanking_3: str = ""
    specificity_value: float = 0.0


@dataclass
class HomologGroup:
    """A/B/D 同源序列组"""
    target_genome: str
    target_seq: str
    target_accession: str
    homologs: Dict[str, str] = field(default_factory=dict)
    homolog_accessions: Dict[str, str] = field(default_factory=dict)
    snp_sites: List[SNPSite] = field(default_factory=list)
    alignment_identity: Dict[str, float] = field(default_factory=dict)


@dataclass
class PrimerCandidate:
    """候选引物"""
    sequence: str
    direction: str
    global_position: int
    length: int
    tm: float
    gc: float
    quality_score: float
    region_specificity: float
    blast_specificity: Optional[float] = None
    blast_result: Optional[BlastResult] = None
    three_prime_snp: bool = False
    snp_count_in_primer: int = 0
    final_score: float = 0.0


# ============================================================
# BLAST 结果缓存
# ============================================================

_BLAST_CACHE: Dict[str, BlastResult] = {}


def _cache_key(sequence: str, entrez_query: str, database: str) -> str:
    raw = f"{sequence.upper()}|{entrez_query}|{database}"
    return hashlib.md5(raw.encode()).hexdigest()


def clear_blast_cache():
    _BLAST_CACHE.clear()


# ============================================================
# 亚基因组识别（增强版）
# ============================================================

def identify_wheat_subgenome(description: str) -> Tuple[str, str]:
    """从 BLAST 命中描述中识别小麦亚基因组 (A/B/D)"""
    desc = description.upper()

    m = re.search(r'TRAES\w{0,6}?(\d)([ABD])\d{2}', desc)
    if m:
        return m.group(2), f"Chr{m.group(1)}{m.group(2)}"

    m = re.search(r'\bCHR(\d)([ABD])(?:_|\b)', desc)
    if m:
        return m.group(2), f"Chr{m.group(1)}{m.group(2)}"

    m = re.search(r'CHROMOSOME\s*(\d)\s*([ABD])\b', desc)
    if m:
        return m.group(2), f"Chr{m.group(1)}{m.group(2)}"

    m = re.search(r'\b(\d)([ABD])\b', desc)
    if m:
        return m.group(2), f"Chr{m.group(1)}{m.group(2)}"

    for kw, genome in [('GENOME A', 'A'), ('GENOME B', 'B'), ('GENOME D', 'D'),
                        ('GENOME_A', 'A'), ('GENOME_B', 'B'), ('GENOME_D', 'D'),
                        ('SUBGENOME A', 'A'), ('SUBGENOME B', 'B'), ('SUBGENOME D', 'D')]:
        if kw in desc:
            return genome, 'Unknown'

    if 'AEGILOPS TAUSCHII' in desc:
        return 'D', 'Donor'
    if 'TRITICUM URARTU' in desc:
        return 'A', 'Donor'

    return 'Unknown', 'Unknown'


# ============================================================
# 解析 BLAST XML 记录（增强版）
# ============================================================

def _parse_blast_record(record, query_length: int, max_hits: int = 50,
                        query_seq: str = "") -> BlastResult:
    """将 Biopython BLAST XML record 解析为 BlastResult"""
    query_label = (record.query or "unknown")[:80]
    hits: List[BlastHit] = []
    subgenome_counter: Dict[str, int] = {'A': 0, 'B': 0, 'D': 0, 'Unknown': 0}

    for alignment in record.alignments[:max_hits]:
        if not alignment.hsps:
            continue

        best_hsp = max(alignment.hsps, key=lambda h: h.score)
        identity_pct = best_hsp.identities / max(best_hsp.align_length, 1) * 100
        coverage_pct = (abs(best_hsp.query_end - best_hsp.query_start) + 1) / max(query_length, 1) * 100

        subgenome, chromosome = identify_wheat_subgenome(alignment.title)
        subgenome_counter[subgenome] = subgenome_counter.get(subgenome, 0) + 1

        hsps_info: List[HSPInfo] = []
        all_mismatch_positions: List[int] = []
        three_prime_ok = True

        for hsp in alignment.hsps:
            q_aligned = str(hsp.query) if hasattr(hsp, 'query') else ""
            s_aligned = str(hsp.sbjct) if hasattr(hsp, 'sbjct') else ""
            mid = str(hsp.match) if hasattr(hsp, 'match') else ""

            if q_aligned and s_aligned:
                for k, (qb, sb) in enumerate(zip(q_aligned, s_aligned)):
                    if qb != sb and qb != '-' and sb != '-':
                        abs_pos = (hsp.query_start - 1) + k
                        all_mismatch_positions.append(abs_pos)

                if hsp.query_end >= query_length - 2:
                    tail = q_aligned[-5:]
                    stail = s_aligned[-5:]
                    if tail != stail:
                        three_prime_ok = False

            hsps_info.append(HSPInfo(
                query_start=hsp.query_start - 1,
                query_end=hsp.query_end,
                subject_start=hsp.sbjct_start - 1,
                subject_end=hsp.sbjct_end,
                identity=hsp.identities / max(hsp.align_length, 1) * 100,
                identity_count=hsp.identities,
                align_length=hsp.align_length,
                mismatches=hsp.align_length - hsp.identities,
                gaps=hsp.gaps if hasattr(hsp, 'gaps') and hsp.gaps else 0,
                score=float(hsp.score),
                e_value=float(hsp.expect),
                query_seq=q_aligned,
                subject_seq=s_aligned,
                midline=mid,
            ))

        hits.append(BlastHit(
            title=alignment.title[:200],
            accession=getattr(alignment, 'accession', ''),
            subject_length=alignment.length,
            best_score=float(best_hsp.score),
            best_e_value=float(best_hsp.expect),
            best_identity=round(identity_pct, 1),
            query_coverage=round(min(coverage_pct, 100.0), 1),
            subgenome=subgenome,
            chromosome=chromosome,
            hsps=hsps_info,
            three_prime_match=three_prime_ok,
            mismatch_positions=sorted(set(all_mismatch_positions)),
        ))

    specificity_score, is_specific, warnings, suggestions = _evaluate_specificity(
        hits, query_length
    )

    return BlastResult(
        query=query_label,
        query_length=query_length,
        database="NCBI nt",
        total_hits=len(record.alignments),
        hits=hits,
        is_specific=is_specific,
        specificity_score=round(specificity_score, 1),
        subgenome_hits=subgenome_counter,
        warnings=warnings,
        suggestions=suggestions,
    )


def _evaluate_specificity(hits: List[BlastHit], query_length: int,
                           ) -> Tuple[float, bool, List[str], List[str]]:
    """综合评估特异性"""
    warnings: List[str] = []
    suggestions: List[str] = []

    perfect_hits = [h for h in hits if h.best_identity >= 95 and h.query_coverage >= 70]
    high_hits = [h for h in hits if h.best_identity >= 85 and h.query_coverage >= 50]

    if not hits:
        return 0.0, False, ["未找到任何命中，请检查序列或网络连接"], ["尝试更宽松的搜索或检查序列是否正确"]

    if not perfect_hits and not high_hits:
        return 0.0, False, [
            "未找到高相似度命中（<85%），序列可能不在小麦数据库中"
        ], ["检查序列物种来源是否正确"]

    genome_perfect: Dict[str, List[BlastHit]] = {'A': [], 'B': [], 'D': [], 'Unknown': []}
    for h in perfect_hits:
        genome_perfect[h.subgenome].append(h)

    is_primer = query_length < 50

    if is_primer:
        return _evaluate_primer_specificity(hits, perfect_hits, high_hits, genome_perfect,
                                             warnings, suggestions)
    else:
        return _evaluate_gene_specificity(hits, perfect_hits, high_hits, genome_perfect,
                                           warnings, suggestions)


def _evaluate_primer_specificity(
    hits, perfect_hits, high_hits, genome_perfect,
    warnings, suggestions
) -> Tuple[float, bool, List[str], List[str]]:
    """引物级别特异性评估"""
    score = 100.0
    is_specific = True

    n_perfect = len(perfect_hits)
    if n_perfect == 0:
        score = 60.0
        warnings.append("未找到完美命中（≥95%），引物可能与目标不匹配")
        suggestions.append("检查引物序列是否正确，或尝试延长引物")
        is_specific = False
    elif n_perfect == 1:
        score = 100.0
        warnings.append("✅ 仅有1个完美命中，特异性优秀")
    elif n_perfect <= 3:
        genomes = {h.subgenome for h in perfect_hits if h.subgenome != 'Unknown'}
        if len(genomes) <= 3 and genomes.issubset({'A', 'B', 'D'}):
            score = 85.0
            is_specific = True
            warnings.append(
                f"在 {'/'.join(sorted(genomes))} 亚基因组共 {n_perfect} 个完美命中（六倍体正常）"
            )
            three_prime_mismatched = [h for h in perfect_hits if not h.three_prime_match]
            if three_prime_mismatched:
                score += 5
                suggestions.append("部分同源命中3'端有错配，有助于特异扩增")
        else:
            score = 55.0
            is_specific = False
            warnings.append(f"同一亚基因组有多个完美命中，可能扩增多个位点")
            suggestions.append("尝试选择3'端落在亚基因组特异SNP上的引物")
    else:
        score = max(0.0, 100.0 - (n_perfect - 1) * 12)
        is_specific = False
        warnings.append(f"存在 {n_perfect} 个完美命中（≥95%），特异性差")
        suggestions.append("重新设计引物，确保3'端落在独特序列区域")

    near_miss = [h for h in high_hits if h not in perfect_hits and h.best_identity >= 90]
    if near_miss:
        score -= min(15, len(near_miss) * 3)
        warnings.append(f"另有 {len(near_miss)} 个 90-95% 相似度命中（潜在非特异扩增风险）")

    return max(0, min(100, score)), is_specific, warnings, suggestions


def _evaluate_gene_specificity(
    hits, perfect_hits, high_hits, genome_perfect,
    warnings, suggestions
) -> Tuple[float, bool, List[str], List[str]]:
    """基因级别特异性评估"""
    score = 80.0
    is_specific = True
    genomes_hit = {h.subgenome for h in perfect_hits if h.subgenome != 'Unknown'}

    if len(genomes_hit) <= 3 and genomes_hit.issubset({'A', 'B', 'D'}):
        score = 85.0
        if len(genomes_hit) == 1:
            score = 95.0
            warnings.append(f"✅ 仅在 {list(genomes_hit)[0]} 亚基因组找到同源序列")
        elif len(genomes_hit) == 3:
            warnings.append("在 A/B/D 三个亚基因组均找到同源序列（典型六倍体基因）")
            suggestions.append("可通过同源比对找到亚基因组间差异位点，用于设计特异引物")
        else:
            warnings.append(f"在 {'/'.join(sorted(genomes_hit))} 亚基因组找到同源序列")

        for g in genomes_hit:
            if len(genome_perfect[g]) > 2:
                score -= 10
                warnings.append(f"{g} 亚基因组有 {len(genome_perfect[g])} 个高相似度命中，可能存在基因家族扩增")
    else:
        unknown_count = len(genome_perfect.get('Unknown', []))
        if unknown_count > 5:
            score = 50.0
            is_specific = False
            warnings.append(f"有 {unknown_count} 个无法归类的命中，可能是高重复序列")
            suggestions.append("该序列可能位于重复区域，不适合直接设计特异引物")

    return max(0, min(100, score)), is_specific, warnings, suggestions


# ============================================================
# 核心 BLAST 函数
# ============================================================

def blast_sequence_ncbi(
    sequence: str,
    entrez_query: str = "Triticum aestivum[Organism]",
    database: str = "nt",
    hitlist: int = 30,
    expect: float = 10.0,
    email: str = "kasp_tool@research.edu",
    megablast: bool = False,
    use_cache: bool = True,
    max_retries: int = 2,
) -> Optional[BlastResult]:
    """NCBI 远程 BLAST（带缓存和重试）"""
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("需要安装 biopython：pip install biopython")

    NCBIWWW, NCBIXML, Entrez = _import_bio()
    Entrez.email = email

    sequence = re.sub(r'[^ATGCatgcNn]', '', sequence)
    sequence = re.sub(r'\s+', '', sequence).upper()

    if len(sequence) < 10:
        raise ValueError(f"序列过短（{len(sequence)}bp），至少需要 10bp")

    ck = _cache_key(sequence, entrez_query, database)
    if use_cache and ck in _BLAST_CACHE:
        return _BLAST_CACHE[ck]

    blast_kwargs = dict(
        program="blastn",
        database=database,
        sequence=sequence,
        entrez_query=entrez_query,
        hitlist_size=hitlist,
        expect=expect,
        format_type="XML",
    )

    if len(sequence) < 50:
        blast_kwargs["word_size"] = 7
    elif megablast:
        blast_kwargs["megablast"] = "on"

    t0 = time.time()
    last_error = None

    for attempt in range(max_retries + 1):
        try:
            result_handle = NCBIWWW.qblast(**blast_kwargs)
            records = list(NCBIXML.parse(result_handle))
            if not records:
                return None
            result = _parse_blast_record(
                records[0], query_length=len(sequence),
                max_hits=hitlist, query_seq=sequence
            )
            result.elapsed_seconds = round(time.time() - t0, 1)
            if use_cache:
                _BLAST_CACHE[ck] = result
            return result
        except Exception as e:
            last_error = e
            if attempt < max_retries:
                time.sleep(5 * (attempt + 1))

    raise RuntimeError(f"BLAST 失败（{max_retries + 1} 次尝试）：{last_error}")


# ============================================================
# 本地 BLAST 支持
# ============================================================

def check_local_blast_available() -> Tuple[bool, str]:
    """检查 BLAST+ 命令行工具是否可用"""
    try:
        r = subprocess.run(["blastn", "-version"], capture_output=True, text=True, timeout=5)
        if r.returncode == 0:
            version = r.stdout.strip().split('\n')[0]
            return True, version
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return False, ""


def build_local_database(fasta_path: str, db_name: str = "wheat_db",
                          output_dir: str = None) -> Tuple[bool, str]:
    """使用 makeblastdb 建立本地数据库"""
    if not os.path.isfile(fasta_path):
        return False, f"文件不存在: {fasta_path}"

    if output_dir is None:
        output_dir = os.path.dirname(fasta_path)
    db_path = os.path.join(output_dir, db_name)

    try:
        r = subprocess.run(
            ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl",
             "-out", db_path, "-title", db_name],
            capture_output=True, text=True, timeout=600,
        )
        if r.returncode == 0:
            return True, db_path
        return False, r.stderr
    except FileNotFoundError:
        return False, "makeblastdb 未找到，请安装 NCBI BLAST+"
    except subprocess.TimeoutExpired:
        return False, "建库超时（>10 分钟）"


def blast_local(
    sequence: str, db_path: str, hitlist: int = 30, evalue: float = 10.0
) -> Optional[BlastResult]:
    """使用本地 BLAST 数据库比对"""
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("需要安装 biopython")

    _, NCBIXML, _ = _import_bio()

    sequence = re.sub(r'[^ATGCatgcNn]', '', sequence).upper()
    if len(sequence) < 10:
        raise ValueError("序列过短")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">query\n{sequence}\n")
        query_file = f.name

    out_file = query_file + ".xml"

    try:
        word_size = 7 if len(sequence) < 50 else 11
        r = subprocess.run(
            ["blastn", "-query", query_file, "-db", db_path,
             "-outfmt", "5", "-out", out_file,
             "-max_target_seqs", str(hitlist),
             "-evalue", str(evalue),
             "-word_size", str(word_size)],
            capture_output=True, text=True, timeout=120,
        )
        if r.returncode != 0:
            return None

        with open(out_file) as fh:
            records = list(NCBIXML.parse(fh))
        if not records:
            return None
        result = _parse_blast_record(records[0], len(sequence), hitlist, sequence)
        result.database = f"Local: {os.path.basename(db_path)}"
        return result
    finally:
        for p in (query_file, out_file):
            try:
                os.unlink(p)
            except OSError:
                pass


# ============================================================
# 同源序列比对 + 差异位点检测
# ============================================================

def extract_homolog_sequences(blast_result: BlastResult, query_seq: str,
                               target_genome: str = "A") -> Optional[HomologGroup]:
    """
    从 BLAST 结果中提取 A/B/D 同源序列，进行 Pairwise 比对，找出差异位点。
    """
    if not blast_result.hits:
        return None

    best_per_genome: Dict[str, BlastHit] = {}
    for hit in blast_result.hits:
        g = hit.subgenome
        if g in ('A', 'B', 'D'):
            if g not in best_per_genome or hit.best_score > best_per_genome[g].best_score:
                best_per_genome[g] = hit

    if target_genome not in best_per_genome:
        return None

    target_hit = best_per_genome[target_genome]

    group = HomologGroup(
        target_genome=target_genome,
        target_seq=query_seq,
        target_accession=target_hit.accession,
    )

    for genome, hit in best_per_genome.items():
        if genome == target_genome:
            continue
        best_hsp = max(hit.hsps, key=lambda h: h.score)
        if not best_hsp.query_seq or not best_hsp.subject_seq:
            continue

        group.homologs[genome] = best_hsp.subject_seq
        group.homolog_accessions[genome] = hit.accession
        group.alignment_identity[genome] = best_hsp.identity

    if not group.homologs:
        return group

    # 找差异位点
    snp_sites: List[SNPSite] = []

    for genome, hit in best_per_genome.items():
        if genome == target_genome:
            continue
        best_hsp = max(hit.hsps, key=lambda h: h.score)
        q_seq = best_hsp.query_seq
        s_seq = best_hsp.subject_seq

        if not q_seq or not s_seq:
            continue

        q_pos = best_hsp.query_start
        for i, (qb, sb) in enumerate(zip(q_seq, s_seq)):
            if qb == '-':
                continue
            if sb == '-':
                snp = SNPSite(
                    position=q_pos, ref_base=qb,
                    alt_bases={genome: '-'}, is_indel=True,
                )
                snp_sites.append(snp)
            elif qb.upper() != sb.upper():
                snp = SNPSite(
                    position=q_pos, ref_base=qb.upper(),
                    alt_bases={genome: sb.upper()},
                )
                snp_sites.append(snp)
            q_pos += 1

    # 合并同位置差异
    merged: Dict[int, SNPSite] = {}
    for snp in snp_sites:
        if snp.position in merged:
            merged[snp.position].alt_bases.update(snp.alt_bases)
            if snp.is_indel:
                merged[snp.position].is_indel = True
        else:
            snp.flanking_5 = query_seq[max(0, snp.position - 10): snp.position]
            snp.flanking_3 = query_seq[snp.position + 1: snp.position + 11]
            merged[snp.position] = snp

    # 评估 SNP 价值
    for pos, snp in merged.items():
        value = 50.0
        if len(snp.alt_bases) >= 2:
            value += 20
        if not snp.is_indel:
            value += 10
        flanking = snp.flanking_5 + snp.flanking_3
        if not re.search(r'(.)\1{3,}', flanking):
            value += 10
        gc = (flanking.count('G') + flanking.count('C')) / max(len(flanking), 1) * 100
        if 30 <= gc <= 65:
            value += 10
        snp.specificity_value = min(100, value)

    group.snp_sites = sorted(merged.values(), key=lambda s: s.specificity_value, reverse=True)
    return group


# ============================================================
# 特异性区域图谱（v2 增强）
# ============================================================

def analyze_specificity_for_design(
    gene_sequence: str,
    blast_result: BlastResult,
    target_genome: str = "A",
    min_region_len: int = 20,
    homolog_group: Optional[HomologGroup] = None,
) -> List[Dict]:
    """
    构建 per-position 特异性图谱。
    v2：整合 SNP 位点信息，包含 SNP 的区域排名更高。
    """
    gene_len = len(gene_sequence)

    target_cov = [0.0] * gene_len
    off_cov = [0.0] * gene_len

    for hit in blast_result.hits:
        for hsp in hit.hsps:
            s = max(0, hsp.query_start)
            e = min(gene_len, hsp.query_end)
            ident = hsp.identity
            if hit.subgenome == target_genome:
                for pos in range(s, e):
                    target_cov[pos] = max(target_cov[pos], ident)
            elif hit.subgenome in ('A', 'B', 'D'):
                for pos in range(s, e):
                    off_cov[pos] = max(off_cov[pos], ident)

    # SNP 加分
    snp_bonus = [0.0] * gene_len
    if homolog_group and homolog_group.snp_sites:
        for snp in homolog_group.snp_sites:
            p = snp.position
            if 0 <= p < gene_len:
                bonus = 15.0 if len(snp.alt_bases) >= 2 else 10.0
                for offset in range(-5, 6):
                    idx = p + offset
                    if 0 <= idx < gene_len:
                        dist_factor = 1.0 - abs(offset) / 6.0
                        snp_bonus[idx] = max(snp_bonus[idx], bonus * dist_factor)

    spec = []
    for i in range(gene_len):
        tgt = target_cov[i]
        off = off_cov[i]
        if tgt >= 90:
            if off <= 70:
                base_score = 95.0
            elif off <= 85:
                base_score = 80.0 - (off - 70) * 1.5
            elif off <= 95:
                base_score = 55.0 - (off - 85) * 2.0
            else:
                base_score = 30.0
        elif tgt >= 70:
            base_score = 60.0
        else:
            base_score = 40.0
        spec.append(round(max(0.0, min(100.0, base_score + snp_bonus[i])), 1))

    regions = []
    i = 0
    threshold = 65.0
    while i < gene_len:
        if spec[i] >= threshold:
            j = i + 1
            while j < gene_len and spec[j] >= threshold:
                j += 1
            length = j - i
            if length >= min_region_len:
                avg_score = sum(spec[i:j]) / length
                t_max = max(target_cov[i:j])
                o_max = max(off_cov[i:j]) if any(off_cov[i:j]) else 0
                snp_max = max(snp_bonus[i:j]) if any(snp_bonus[i:j]) else 0

                snp_count = 0
                if homolog_group:
                    for snp in homolog_group.snp_sites:
                        if i <= snp.position < j:
                            snp_count += 1

                regions.append({
                    'start': i,
                    'end': j,
                    'length': length,
                    'score': round(avg_score, 1),
                    'target_coverage': round(t_max, 1),
                    'off_target_coverage': round(o_max, 1),
                    'snp_count': snp_count,
                    'has_snp_bonus': snp_max > 0,
                    'reason': (
                        f"目标({target_genome})覆盖 {t_max:.0f}%，"
                        f"非目标覆盖 {o_max:.0f}%，"
                        f"含 {snp_count} 个 SNP，"
                        f"特异性分 {avg_score:.0f}"
                    ),
                    'sequence': gene_sequence[i:j],
                    'specificity_profile': spec[i:j],
                })
            i = j
        else:
            i += 1

    regions.sort(key=lambda r: (r['snp_count'] > 0, r['score']), reverse=True)
    return regions
