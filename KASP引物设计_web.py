"""
KASP & 常规PCR 引物设计工具 - Streamlit Web版
版本: v5.0 Web
功能: KASP引物设计、常规PCR引物设计、质量评估、CSV导出
"""

import streamlit as st
import re
import csv
import io
from datetime import datetime
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

# ==================== 页面配置 ====================
st.set_page_config(
    page_title="引物设计工具 v5.0",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ==================== 自定义样式 ====================
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1E88E5;
        text-align: center;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .primer-box {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
        font-family: monospace;
    }
    .score-excellent { color: #4CAF50; font-weight: bold; }
    .score-good { color: #8BC34A; font-weight: bold; }
    .score-ok { color: #FFC107; font-weight: bold; }
    .score-poor { color: #FF9800; font-weight: bold; }
    .score-bad { color: #f44336; font-weight: bold; }
    .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 1rem 0;
    }
    .success-box {
        background-color: #d4edda;
        border: 1px solid #28a745;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 1rem 0;
    }
    .info-box {
        background-color: #e7f3ff;
        border: 1px solid #1E88E5;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

# ==================== 配置类 ====================
@dataclass
class KASPConfig:
    """KASP引物设计配置参数"""
    FAM_TAIL: str = "GAAGGTGACCAAGTTCATGCT"
    HEX_TAIL: str = "GAAGGTCGGAGTCAACGGATT"
    MIN_PRIMER_LEN: int = 18
    MAX_PRIMER_LEN: int = 30
    OPTIMAL_PRIMER_LEN: int = 20
    MIN_TM: float = 55.0
    MAX_TM: float = 68.0
    OPTIMAL_TM: float = 62.0
    MIN_GC: float = 35.0
    MAX_GC: float = 65.0
    OPTIMAL_GC_MIN: float = 40.0
    OPTIMAL_GC_MAX: float = 60.0
    MAX_TM_DIFF: float = 2.0
    REV_MIN_DISTANCE: int = 50
    REV_MAX_DISTANCE: int = 150
    PRODUCT_MIN: int = 80
    PRODUCT_MAX: int = 200
    MISMATCH_POSITIONS: List[int] = None
    
    def __post_init__(self):
        if self.MISMATCH_POSITIONS is None:
            self.MISMATCH_POSITIONS = [-2, -3, -4]


@dataclass
class RegularPCRConfig:
    """常规PCR引物设计配置"""
    MIN_PRIMER_LEN: int = 18
    MAX_PRIMER_LEN: int = 25
    OPTIMAL_PRIMER_LEN: int = 20
    MIN_TM: float = 55.0
    MAX_TM: float = 68.0
    OPTIMAL_TM: float = 60.0
    MIN_GC: float = 40.0
    MAX_GC: float = 60.0
    OPTIMAL_GC_MIN: float = 45.0
    OPTIMAL_GC_MAX: float = 55.0
    MAX_TM_DIFF: float = 2.0
    PRODUCT_MIN: int = 150
    PRODUCT_MAX: int = 500


# ==================== 核心计算函数 ====================

def calc_gc_content(seq: str) -> float:
    """计算GC含量百分比"""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def calc_tm_nearest_neighbor(seq: str, na_conc: float = 50.0) -> float:
    """
    使用简化的经验公式计算Tm值
    """
    seq = seq.upper()
    length = len(seq)
    
    if length < 14:
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        tm = 2 * at_count + 4 * gc_count
    else:
        gc_percent = calc_gc_content(seq)
        tm = 81.5 + 0.41 * gc_percent - 675 / length
        salt_correction = 16.6 * (0.69897 + (-3.0))  # log10(0.001) ≈ -3
        tm = tm + salt_correction * 0.1
    
    return round(tm, 1)


def reverse_complement(seq: str) -> str:
    """生成反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def check_hairpin(seq: str, min_stem: int = 4, min_loop: int = 3) -> bool:
    """检测发夹结构"""
    seq = seq.upper()
    length = len(seq)
    
    for i in range(length - min_stem - min_loop):
        for stem_len in range(min_stem, min(8, (length - i - min_loop) // 2 + 1)):
            stem1 = seq[i:i + stem_len]
            for loop_len in range(min_loop, min(10, length - i - 2 * stem_len + 1)):
                j = i + stem_len + loop_len
                if j + stem_len <= length:
                    stem2 = seq[j:j + stem_len]
                    if stem1 == reverse_complement(stem2)[::-1]:
                        return True
    return False


def check_self_dimer(seq: str, min_complementary: int = 4) -> bool:
    """检测自身二聚体"""
    seq = seq.upper()
    rc = reverse_complement(seq)
    length = len(seq)
    
    for i in range(length - min_complementary + 1):
        segment = seq[i:i + min_complementary]
        if segment in rc:
            return True
    
    end_seq = seq[-6:]
    if any(end_seq[i:i+3] in rc[-10:] for i in range(4)):
        return True
    
    return False


def check_primer_dimer(seq1: str, seq2: str, min_complementary: int = 4) -> bool:
    """检测引物二聚体"""
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    rc2 = reverse_complement(seq2)
    
    for i in range(len(seq1) - min_complementary + 1):
        segment = seq1[i:i + min_complementary]
        if segment in rc2:
            return True
    
    end1 = seq1[-6:]
    end2_rc = reverse_complement(seq2[-6:])
    for i in range(4):
        if end1[i:i+3] in end2_rc:
            return True
    
    return False


def check_3prime_stability(seq: str) -> Tuple[bool, str]:
    """检查3'端稳定性"""
    end_5 = seq[-5:].upper()
    gc_count = end_5.count('G') + end_5.count('C')
    
    if gc_count > 3:
        return False, "3'端GC过多，可能导致非特异性结合"
    elif gc_count < 1:
        return False, "3'端GC过少，结合不稳定"
    
    if seq[-1].upper() in ['G', 'C']:
        return True, "3'端以G/C结尾，良好"
    else:
        return True, "3'端以A/T结尾，可接受"


def get_strong_mismatch(original_base: str) -> str:
    """获取强错配碱基"""
    strong_mismatches = {
        'A': 'A', 'T': 'T', 'G': 'A', 'C': 'A'
    }
    return strong_mismatches.get(original_base.upper(), 'A')


def evaluate_primer_quality(seq: str, config=None) -> Dict:
    """综合评估引物质量"""
    if config is None:
        config = KASPConfig()
    
    seq = seq.upper()
    result = {
        'sequence': seq,
        'length': len(seq),
        'gc_content': calc_gc_content(seq),
        'tm': calc_tm_nearest_neighbor(seq),
        'has_hairpin': check_hairpin(seq),
        'has_self_dimer': check_self_dimer(seq),
        'three_prime_ok': check_3prime_stability(seq)[0],
        'three_prime_msg': check_3prime_stability(seq)[1],
        'issues': [],
        'score': 100
    }
    
    # 评分扣分逻辑
    if result['tm'] < config.MIN_TM:
        result['issues'].append(f"Tm过低({result['tm']}°C)")
        result['score'] -= 15
    elif result['tm'] > config.MAX_TM:
        result['issues'].append(f"Tm过高({result['tm']}°C)")
        result['score'] -= 10
    elif abs(result['tm'] - config.OPTIMAL_TM) > 3:
        result['issues'].append(f"Tm偏离最优值({result['tm']}°C)")
        result['score'] -= 5
    
    if result['gc_content'] < config.MIN_GC:
        result['issues'].append(f"GC含量过低({result['gc_content']:.1f}%)")
        result['score'] -= 15
    elif result['gc_content'] > config.MAX_GC:
        result['issues'].append(f"GC含量过高({result['gc_content']:.1f}%)")
        result['score'] -= 15
    elif not (config.OPTIMAL_GC_MIN <= result['gc_content'] <= config.OPTIMAL_GC_MAX):
        result['score'] -= 5
    
    if result['has_hairpin']:
        result['issues'].append("可能形成发夹结构")
        result['score'] -= 10
    
    if result['has_self_dimer']:
        result['issues'].append("自身二聚体风险")
        result['score'] -= 10
    
    if not result['three_prime_ok']:
        result['issues'].append(result['three_prime_msg'])
        result['score'] -= 10
    
    if len(seq) < config.MIN_PRIMER_LEN or len(seq) > config.MAX_PRIMER_LEN:
        result['issues'].append(f"长度不在最优范围({len(seq)}bp)")
        result['score'] -= 5
    
    result['score'] = max(0, result['score'])
    
    return result


def get_quality_grade(score: float) -> Tuple[str, str, str]:
    """根据评分返回质量等级、星级和CSS类"""
    if score >= 85:
        return "优秀", "★★★★★", "score-excellent"
    elif score >= 75:
        return "良好", "★★★★☆", "score-good"
    elif score >= 65:
        return "合格", "★★★☆☆", "score-ok"
    elif score >= 55:
        return "可用", "★★☆☆☆", "score-poor"
    elif score >= 45:
        return "勉强", "★☆☆☆☆", "score-bad"
    else:
        return "需优化", "☆☆☆☆☆", "score-bad"


# ==================== KASP引物设计 ====================

def parse_snp_sequence(seq_with_snp: str) -> Tuple[str, str, str, str]:
    """解析包含SNP标记的序列"""
    seq_with_snp = seq_with_snp.strip().upper()
    seq_with_snp = re.sub(r'\s+', '', seq_with_snp)
    
    pattern = r'\[([ATGC])/([ATGC])\]'
    match = re.search(pattern, seq_with_snp, re.IGNORECASE)
    
    if not match:
        raise ValueError("未找到SNP标记，请使用格式如 [A/G] 或 [C/T]")
    
    allele1 = match.group(1).upper()
    allele2 = match.group(2).upper()
    
    snp_pos = match.start()
    upstream = seq_with_snp[:snp_pos]
    downstream = seq_with_snp[match.end():]
    
    upstream = re.sub(r'[^ATGC]', '', upstream.upper())
    downstream = re.sub(r'[^ATGC]', '', downstream.upper())
    
    return upstream, downstream, allele1, allele2


def design_kasp_primers_multi(upstream: str, downstream: str, allele1: str, allele2: str, 
                              config: KASPConfig = None, num_schemes: int = 5) -> List[Dict]:
    """设计多套KASP引物方案"""
    if config is None:
        config = KASPConfig()
    
    all_schemes = []
    
    # 生成不同长度的正向引物
    for primer_len in range(config.MIN_PRIMER_LEN, min(config.MAX_PRIMER_LEN + 1, len(upstream) + 1)):
        core_seq = upstream[-(primer_len - 1):]
        
        if len(core_seq) < config.MIN_PRIMER_LEN - 1:
            continue
        
        for mismatch_pos in config.MISMATCH_POSITIONS:
            if abs(mismatch_pos) >= len(core_seq):
                continue
            
            mismatch_idx = len(core_seq) + mismatch_pos
            original_base = core_seq[mismatch_idx]
            mismatch_base = get_strong_mismatch(original_base)
            
            if mismatch_base == original_base:
                continue
            
            core_with_mismatch = core_seq[:mismatch_idx] + mismatch_base + core_seq[mismatch_idx + 1:]
            
            fwd_allele1 = core_with_mismatch + allele1
            fwd_allele2 = core_with_mismatch + allele2
            
            fwd_with_fam = config.FAM_TAIL + fwd_allele1
            fwd_with_hex = config.HEX_TAIL + fwd_allele2
            
            # 评估正向引物
            eval1 = evaluate_primer_quality(fwd_allele1, config)
            eval2 = evaluate_primer_quality(fwd_allele2, config)
            
            tm_diff = abs(eval1['tm'] - eval2['tm'])
            
            # 搜索反向引物
            for rev_dist in range(config.REV_MIN_DISTANCE, min(config.REV_MAX_DISTANCE + 1, len(downstream) - config.MIN_PRIMER_LEN + 1)):
                for rev_len in range(config.MIN_PRIMER_LEN, min(config.MAX_PRIMER_LEN + 1, len(downstream) - rev_dist + 1)):
                    rev_start = rev_dist
                    rev_end = rev_dist + rev_len
                    
                    if rev_end > len(downstream):
                        continue
                    
                    rev_seq = reverse_complement(downstream[rev_start:rev_end])
                    eval_rev = evaluate_primer_quality(rev_seq, config)
                    
                    # 检查引物二聚体
                    has_dimer = (check_primer_dimer(fwd_allele1, rev_seq) or 
                                check_primer_dimer(fwd_allele2, rev_seq))
                    
                    # 计算综合评分
                    avg_fwd_score = (eval1['score'] + eval2['score']) / 2
                    total_score = (avg_fwd_score * 0.5 + eval_rev['score'] * 0.3)
                    
                    if tm_diff <= 1.0:
                        total_score += 15
                    elif tm_diff <= 2.0:
                        total_score += 5
                    else:
                        total_score -= 10
                    
                    if has_dimer:
                        total_score -= 10
                    
                    product_size = len(upstream) + 1 + rev_dist + rev_len
                    if config.PRODUCT_MIN <= product_size <= config.PRODUCT_MAX:
                        total_score += 5
                    
                    total_score = max(0, min(100, total_score))
                    
                    scheme = {
                        'fwd_allele1_full': fwd_with_fam,
                        'fwd_allele2_full': fwd_with_hex,
                        'fwd_allele1_core': fwd_allele1,
                        'fwd_allele2_core': fwd_allele2,
                        'reverse': rev_seq,
                        'allele1': allele1,
                        'allele2': allele2,
                        'mismatch_pos': mismatch_pos,
                        'mismatch_change': f"{original_base}→{mismatch_base}",
                        'eval_fwd1': eval1,
                        'eval_fwd2': eval2,
                        'eval_rev': eval_rev,
                        'tm_diff': tm_diff,
                        'has_dimer': has_dimer,
                        'product_size': product_size,
                        'rev_distance': rev_dist,
                        'total_score': total_score
                    }
                    all_schemes.append(scheme)
    
    # 按评分排序并去重
    all_schemes.sort(key=lambda x: x['total_score'], reverse=True)
    
    unique_schemes = []
    seen = set()
    for scheme in all_schemes:
        key = (scheme['fwd_allele1_core'], scheme['reverse'])
        if key not in seen:
            seen.add(key)
            unique_schemes.append(scheme)
            if len(unique_schemes) >= num_schemes:
                break
    
    return unique_schemes


# ==================== 常规PCR引物设计 ====================

def design_regular_primers(sequence: str, config: RegularPCRConfig = None, 
                          num_pairs: int = 5, target_start: int = None, 
                          target_end: int = None) -> List[Dict]:
    """设计常规PCR引物对"""
    if config is None:
        config = RegularPCRConfig()
    
    sequence = re.sub(r'[^ATGC]', '', sequence.upper())
    seq_len = len(sequence)
    
    if target_start is None:
        target_start = 0
    if target_end is None:
        target_end = seq_len
    
    all_pairs = []
    
    # 搜索正向引物区域
    fwd_search_start = max(0, target_start - 100)
    fwd_search_end = min(target_start + 50, seq_len - config.PRODUCT_MIN)
    
    # 搜索反向引物区域
    rev_search_start = max(target_end - 50, config.PRODUCT_MIN)
    rev_search_end = min(target_end + 100, seq_len)
    
    for fwd_start in range(fwd_search_start, fwd_search_end, 5):
        for fwd_len in range(config.MIN_PRIMER_LEN, config.MAX_PRIMER_LEN + 1):
            fwd_end = fwd_start + fwd_len
            if fwd_end > seq_len:
                continue
            
            fwd_seq = sequence[fwd_start:fwd_end]
            fwd_eval = evaluate_primer_quality(fwd_seq, config)
            
            if fwd_eval['score'] < 40:
                continue
            
            for rev_end in range(rev_search_start, rev_search_end, 5):
                for rev_len in range(config.MIN_PRIMER_LEN, config.MAX_PRIMER_LEN + 1):
                    rev_start = rev_end - rev_len
                    if rev_start < 0:
                        continue
                    
                    product_size = rev_end - fwd_start
                    if not (config.PRODUCT_MIN <= product_size <= config.PRODUCT_MAX):
                        continue
                    
                    rev_seq = reverse_complement(sequence[rev_start:rev_end])
                    rev_eval = evaluate_primer_quality(rev_seq, config)
                    
                    if rev_eval['score'] < 40:
                        continue
                    
                    tm_diff = abs(fwd_eval['tm'] - rev_eval['tm'])
                    has_dimer = check_primer_dimer(fwd_seq, rev_seq)
                    
                    # 计算综合评分
                    total_score = (fwd_eval['score'] + rev_eval['score']) / 2
                    
                    if tm_diff <= 1.0:
                        total_score += 10
                    elif tm_diff <= 2.0:
                        total_score += 5
                    else:
                        total_score -= 10
                    
                    if has_dimer:
                        total_score -= 15
                    
                    if 200 <= product_size <= 400:
                        total_score += 5
                    
                    total_score = max(0, min(100, total_score))
                    
                    pair = {
                        'forward': fwd_seq,
                        'reverse': rev_seq,
                        'fwd_start': fwd_start + 1,
                        'fwd_end': fwd_end,
                        'rev_start': rev_start + 1,
                        'rev_end': rev_end,
                        'fwd_eval': fwd_eval,
                        'rev_eval': rev_eval,
                        'tm_diff': tm_diff,
                        'has_dimer': has_dimer,
                        'product_size': product_size,
                        'total_score': total_score
                    }
                    all_pairs.append(pair)
    
    # 排序并去重
    all_pairs.sort(key=lambda x: x['total_score'], reverse=True)
    
    unique_pairs = []
    seen = set()
    for pair in all_pairs:
        key = (pair['forward'], pair['reverse'])
        if key not in seen:
            seen.add(key)
            unique_pairs.append(pair)
            if len(unique_pairs) >= num_pairs:
                break
    
    return unique_pairs


# ==================== CSV导出函数 ====================

def generate_kasp_csv(schemes: List[Dict], seq_id: str) -> str:
    """生成KASP引物CSV内容"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    writer.writerow(['KASP引物设计报告'])
    writer.writerow(['序列ID', seq_id])
    writer.writerow(['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
    writer.writerow([])
    
    writer.writerow(['方案', '评分', '等级', 
                     'FAM引物(完整)', 'HEX引物(完整)', '通用反向引物',
                     'Allele1-Tm', 'Allele2-Tm', 'Tm差异',
                     'Allele1-GC%', 'Allele2-GC%', 'Rev-Tm', 'Rev-GC%',
                     '产物大小', '错配位置', '错配变化'])
    
    for i, scheme in enumerate(schemes, 1):
        grade, stars, _ = get_quality_grade(scheme['total_score'])
        writer.writerow([
            f"方案{i}", f"{scheme['total_score']:.1f}", f"{grade} {stars}",
            scheme['fwd_allele1_full'], scheme['fwd_allele2_full'], scheme['reverse'],
            f"{scheme['eval_fwd1']['tm']}°C", f"{scheme['eval_fwd2']['tm']}°C", f"{scheme['tm_diff']:.1f}°C",
            f"{scheme['eval_fwd1']['gc_content']:.1f}%", f"{scheme['eval_fwd2']['gc_content']:.1f}%",
            f"{scheme['eval_rev']['tm']}°C", f"{scheme['eval_rev']['gc_content']:.1f}%",
            f"{scheme['product_size']}bp", f"n{scheme['mismatch_pos']}", scheme['mismatch_change']
        ])
    
    return output.getvalue()


def generate_regular_csv(pairs: List[Dict], seq_id: str) -> str:
    """生成常规PCR引物CSV内容"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    writer.writerow(['常规PCR引物设计报告'])
    writer.writerow(['序列ID', seq_id])
    writer.writerow(['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
    writer.writerow([])
    
    writer.writerow(['引物对', '评分', '等级',
                     '正向引物', '反向引物',
                     'Fwd-Tm', 'Rev-Tm', 'Tm差异',
                     'Fwd-GC%', 'Rev-GC%',
                     'Fwd位置', 'Rev位置', '产物大小'])
    
    for i, pair in enumerate(pairs, 1):
        grade, stars, _ = get_quality_grade(pair['total_score'])
        writer.writerow([
            f"引物对{i}", f"{pair['total_score']:.1f}", f"{grade} {stars}",
            pair['forward'], pair['reverse'],
            f"{pair['fwd_eval']['tm']}°C", f"{pair['rev_eval']['tm']}°C", f"{pair['tm_diff']:.1f}°C",
            f"{pair['fwd_eval']['gc_content']:.1f}%", f"{pair['rev_eval']['gc_content']:.1f}%",
            f"{pair['fwd_start']}-{pair['fwd_end']}", f"{pair['rev_start']}-{pair['rev_end']}",
            f"{pair['product_size']}bp"
        ])
    
    return output.getvalue()


# ==================== Streamlit 界面 ====================

def show_kasp_design():
    """KASP引物设计页面"""
    st.markdown("### 🔬 KASP引物设计")
    
    st.markdown("""
    <div class="info-box">
    <b>使用说明：</b><br>
    1. 在输入框中粘贴包含SNP位点的序列<br>
    2. SNP位点使用方括号标记，如 <code>[A/G]</code> 或 <code>[C/T]</code><br>
    3. SNP上游序列建议 ≥50bp，下游序列建议 ≥200bp<br>
    4. 点击"设计引物"按钮开始设计
    </div>
    """, unsafe_allow_html=True)
    
    # 示例序列
    example_seq = """CGTTGCATGAATCCCACAACGCACAGCCGTTGCTCGTCGCCGCCGCCGCCATGGCATTTTTATGTACGCAGAGGAAGAACAAACTCGAG
AAGAAGGCTGAGGAGCTGGAGGAATGGGTCACGGACTA[G/T]GTCGCGATATACTACGCCGACGGACTACTGTCGCGATGGTGATGAA
GGAGACCCGCAAGGCGCTCGGATCGGCTTACCACTCCATGATGATGGTGGAGCAGGTCCACCTGGGGAAGAGCGCCAACTGGGACGAGCTCATCAAC"""
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        seq_input = st.text_area(
            "输入序列（包含SNP标记）",
            value=example_seq,
            height=200,
            help="SNP位点使用 [碱基1/碱基2] 格式标记"
        )
        
        seq_id = st.text_input("序列ID（可选）", value="My_SNP_Marker")
    
    with col2:
        st.markdown("**参数设置**")
        num_schemes = st.slider("生成方案数", 3, 15, 5)
        
        with st.expander("高级参数"):
            min_primer_len = st.number_input("最小引物长度", 15, 25, 18)
            max_primer_len = st.number_input("最大引物长度", 20, 35, 30)
            min_tm = st.number_input("最低Tm (°C)", 50.0, 60.0, 55.0)
            max_tm = st.number_input("最高Tm (°C)", 60.0, 75.0, 68.0)
    
    if st.button("🧬 设计KASP引物", type="primary", use_container_width=True):
        try:
            with st.spinner("正在设计引物..."):
                upstream, downstream, allele1, allele2 = parse_snp_sequence(seq_input)
                
                config = KASPConfig()
                config.MIN_PRIMER_LEN = min_primer_len
                config.MAX_PRIMER_LEN = max_primer_len
                config.MIN_TM = min_tm
                config.MAX_TM = max_tm
                
                schemes = design_kasp_primers_multi(upstream, downstream, allele1, allele2, config, num_schemes)
            
            if not schemes:
                st.error("❌ 未能设计出合适的引物，请检查序列或调整参数")
                return
            
            st.success(f"✅ 成功设计 {len(schemes)} 套引物方案！")
            
            # 显示SNP信息
            st.markdown(f"""
            **SNP信息：** `{allele1}/{allele2}` | 
            **上游序列：** {len(upstream)}bp | 
            **下游序列：** {len(downstream)}bp
            """)
            
            # 显示每个方案
            for i, scheme in enumerate(schemes, 1):
                grade, stars, css_class = get_quality_grade(scheme['total_score'])
                
                with st.expander(f"方案 #{i} | 评分: {scheme['total_score']:.1f} | {grade} {stars}", expanded=(i==1)):
                    col_a, col_b = st.columns(2)
                    
                    with col_a:
                        st.markdown("**正向引物 (等位基因特异性)**")
                        st.code(f"Allele {allele1} (FAM): {scheme['fwd_allele1_full']}")
                        st.caption(f"核心: {scheme['fwd_allele1_core']} | {len(scheme['fwd_allele1_core'])}bp | Tm: {scheme['eval_fwd1']['tm']}°C | GC: {scheme['eval_fwd1']['gc_content']:.1f}%")
                        
                        st.code(f"Allele {allele2} (HEX): {scheme['fwd_allele2_full']}")
                        st.caption(f"核心: {scheme['fwd_allele2_core']} | {len(scheme['fwd_allele2_core'])}bp | Tm: {scheme['eval_fwd2']['tm']}°C | GC: {scheme['eval_fwd2']['gc_content']:.1f}%")
                    
                    with col_b:
                        st.markdown("**反向引物 (通用)**")
                        st.code(f"Common: {scheme['reverse']}")
                        st.caption(f"{len(scheme['reverse'])}bp | Tm: {scheme['eval_rev']['tm']}°C | GC: {scheme['eval_rev']['gc_content']:.1f}% | 距SNP: {scheme['rev_distance']}bp")
                        
                        st.markdown("**产物信息**")
                        tm_status = "✓" if scheme['tm_diff'] <= 1.5 else ("△" if scheme['tm_diff'] <= 2.0 else "✗")
                        st.write(f"Tm差异: {scheme['tm_diff']:.1f}°C {tm_status}")
                        st.write(f"错配位置: n{scheme['mismatch_pos']} ({scheme['mismatch_change']})")
                        st.write(f"产物大小: {scheme['product_size']} bp")
                    
                    # 显示问题
                    all_issues = scheme['eval_fwd1']['issues'] + scheme['eval_fwd2']['issues'] + scheme['eval_rev']['issues']
                    if scheme['has_dimer']:
                        all_issues.append("引物间可能形成二聚体")
                    
                    if all_issues:
                        st.warning("⚠️ 注意事项: " + " | ".join(set(all_issues)))
            
            # 优化建议
            best_score = schemes[0]['total_score'] if schemes else 0
            if best_score < 70:
                st.markdown("---")
                st.markdown("### ⚠️ 质量等级较低 - 优化建议")
                
                suggestions = [
                    ("提供更长的侧翼序列", "上游建议 >100bp，下游建议 >300bp，以便有更多引物设计空间"),
                    ("调整退火温度", "进行梯度PCR优化（52-62°C），找到最佳退火温度"),
                    ("优化反应体系", "使用KASP专用Master Mix，优化Mg²⁺浓度（1.5-2.5mM）"),
                    ("减少引物浓度", "尝试100-200nM引物浓度减少非特异性结合"),
                    ("使用热启动酶", "减少室温下的非特异性扩增")
                ]
                
                for title, detail in suggestions:
                    st.markdown(f"**{title}**")
                    st.caption(detail)
            
            # 导出CSV
            st.markdown("---")
            csv_content = generate_kasp_csv(schemes, seq_id)
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            
            st.download_button(
                label="📥 下载CSV报告",
                data=csv_content,
                file_name=f"{timestamp}-KASP_{seq_id}.csv",
                mime="text/csv",
                use_container_width=True
            )
            
            # 推荐方案快速复制
            st.markdown("---")
            st.markdown("### 📋 推荐方案 (方案 #1)")
            best = schemes[0]
            
            copy_col1, copy_col2, copy_col3 = st.columns(3)
            with copy_col1:
                st.text_input("FAM引物", value=best['fwd_allele1_full'], key="copy_fam")
            with copy_col2:
                st.text_input("HEX引物", value=best['fwd_allele2_full'], key="copy_hex")
            with copy_col3:
                st.text_input("通用引物", value=best['reverse'], key="copy_common")
                
        except ValueError as e:
            st.error(f"❌ 序列解析错误: {e}")
        except Exception as e:
            st.error(f"❌ 设计过程出错: {e}")


def show_regular_pcr_design():
    """常规PCR引物设计页面"""
    st.markdown("### 🧪 常规PCR引物设计")
    
    st.markdown("""
    <div class="info-box">
    <b>使用说明：</b><br>
    1. 输入目标基因/序列（纯碱基序列，无需标记）<br>
    2. 可选择指定目标区域的起止位置<br>
    3. 调整产物大小范围和其他参数<br>
    4. 点击"设计引物"按钮开始设计
    </div>
    """, unsafe_allow_html=True)
    
    example_seq = """ATGGCATTTTTATGTACGCAGAGGAAGAACAAACTCGAGAAGAAGGCTGAGGAGCTGGAGGAATGGGTCACGGACTAGTCGCGATATACTACGCCGACGGACTACTGTCGCGATGGTGATGAAGGAGACCCGCAAGGCGCTCGGATCGGCTTACCACTCCATGATGATGGTGGAGCAGGTCCACCTGGGGAAGAGCGCCAACTGGGACGAGCTCATCAACGAGGTCAAGGCCAAGATCCAGGACAAGGAGGGCATCCCCCCGGACCAGCAGAGGATGATCAACGAGATCAAGATCCTGAACCGCAGGTGA"""
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        seq_input = st.text_area(
            "输入序列",
            value=example_seq,
            height=200,
            help="输入纯碱基序列(A/T/G/C)"
        )
        
        seq_id = st.text_input("序列ID（可选）", value="My_Gene", key="regular_seq_id")
    
    with col2:
        st.markdown("**参数设置**")
        num_pairs = st.slider("生成引物对数", 3, 10, 5, key="regular_num")
        
        st.markdown("**产物大小范围**")
        product_min = st.number_input("最小产物(bp)", 100, 500, 150)
        product_max = st.number_input("最大产物(bp)", 200, 1000, 500)
        
        with st.expander("目标区域（可选）"):
            use_target = st.checkbox("指定目标区域")
            if use_target:
                target_start = st.number_input("起始位置", 1, 10000, 50)
                target_end = st.number_input("结束位置", 1, 10000, 200)
            else:
                target_start = None
                target_end = None
    
    if st.button("🧪 设计PCR引物", type="primary", use_container_width=True):
        try:
            with st.spinner("正在设计引物..."):
                sequence = re.sub(r'[^ATGC]', '', seq_input.upper())
                
                if len(sequence) < product_min + 50:
                    st.error(f"❌ 序列过短（{len(sequence)}bp），请提供更长的序列")
                    return
                
                config = RegularPCRConfig()
                config.PRODUCT_MIN = product_min
                config.PRODUCT_MAX = product_max
                
                pairs = design_regular_primers(
                    sequence, config, num_pairs,
                    target_start if use_target else None,
                    target_end if use_target else None
                )
            
            if not pairs:
                st.error("❌ 未能设计出合适的引物对，请调整参数或检查序列")
                return
            
            st.success(f"✅ 成功设计 {len(pairs)} 对引物！")
            st.caption(f"序列长度: {len(sequence)} bp")
            
            # 显示每对引物
            for i, pair in enumerate(pairs, 1):
                grade, stars, css_class = get_quality_grade(pair['total_score'])
                
                with st.expander(f"引物对 #{i} | 评分: {pair['total_score']:.1f} | {grade} {stars}", expanded=(i==1)):
                    col_a, col_b = st.columns(2)
                    
                    with col_a:
                        st.markdown("**正向引物 (Forward)**")
                        st.code(f"5'- {pair['forward']} -3'")
                        st.caption(f"位置: {pair['fwd_start']}-{pair['fwd_end']} | {len(pair['forward'])}bp | Tm: {pair['fwd_eval']['tm']}°C | GC: {pair['fwd_eval']['gc_content']:.1f}%")
                    
                    with col_b:
                        st.markdown("**反向引物 (Reverse)**")
                        st.code(f"5'- {pair['reverse']} -3'")
                        st.caption(f"位置: {pair['rev_start']}-{pair['rev_end']} | {len(pair['reverse'])}bp | Tm: {pair['rev_eval']['tm']}°C | GC: {pair['rev_eval']['gc_content']:.1f}%")
                    
                    tm_status = "✓" if pair['tm_diff'] <= 1.0 else ("△" if pair['tm_diff'] <= 2.0 else "✗")
                    st.write(f"**Tm差异:** {pair['tm_diff']:.1f}°C {tm_status} | **产物大小:** {pair['product_size']} bp")
                    
                    all_issues = pair['fwd_eval']['issues'] + pair['rev_eval']['issues']
                    if pair['has_dimer']:
                        all_issues.append("引物间可能形成二聚体")
                    
                    if all_issues:
                        st.warning("⚠️ " + " | ".join(set(all_issues)))
            
            # 导出CSV
            st.markdown("---")
            csv_content = generate_regular_csv(pairs, seq_id)
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            
            st.download_button(
                label="📥 下载CSV报告",
                data=csv_content,
                file_name=f"{timestamp}-PCR_{seq_id}.csv",
                mime="text/csv",
                use_container_width=True
            )
            
            # 推荐引物快速复制
            st.markdown("---")
            st.markdown("### 📋 推荐引物对 (#1)")
            best = pairs[0]
            
            copy_col1, copy_col2 = st.columns(2)
            with copy_col1:
                st.text_input("Forward", value=best['forward'], key="copy_fwd")
            with copy_col2:
                st.text_input("Reverse", value=best['reverse'], key="copy_rev")
                
        except Exception as e:
            st.error(f"❌ 设计过程出错: {e}")


def show_primer_analysis():
    """引物分析工具"""
    st.markdown("### 🔍 引物质量分析")
    
    st.markdown("""
    <div class="info-box">
    输入引物序列，快速分析其各项质量指标
    </div>
    """, unsafe_allow_html=True)
    
    primer_input = st.text_input("输入引物序列", placeholder="例如: ATGCGATCGATCGATCG")
    
    if primer_input:
        primer = re.sub(r'[^ATGC]', '', primer_input.upper())
        
        if len(primer) < 10:
            st.warning("引物序列过短")
            return
        
        result = evaluate_primer_quality(primer)
        grade, stars, css_class = get_quality_grade(result['score'])
        
        st.markdown(f"### 分析结果: <span class='{css_class}'>{grade} {stars} ({result['score']:.1f}分)</span>", unsafe_allow_html=True)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("长度", f"{result['length']} bp")
        with col2:
            st.metric("Tm", f"{result['tm']} °C")
        with col3:
            st.metric("GC含量", f"{result['gc_content']:.1f}%")
        with col4:
            st.metric("评分", f"{result['score']}")
        
        st.markdown("---")
        
        check_col1, check_col2 = st.columns(2)
        
        with check_col1:
            st.markdown("**结构检测**")
            st.write(f"{'❌' if result['has_hairpin'] else '✅'} 发夹结构: {'检测到' if result['has_hairpin'] else '未检测到'}")
            st.write(f"{'❌' if result['has_self_dimer'] else '✅'} 自身二聚体: {'有风险' if result['has_self_dimer'] else '无风险'}")
        
        with check_col2:
            st.markdown("**3'端分析**")
            st.write(f"{'✅' if result['three_prime_ok'] else '⚠️'} {result['three_prime_msg']}")
        
        if result['issues']:
            st.warning("⚠️ 发现问题: " + " | ".join(result['issues']))


def show_help():
    """帮助文档"""
    st.markdown("### 📖 使用帮助")
    
    st.markdown("""
    ## KASP引物设计
    
    **KASP (Kompetitive Allele Specific PCR)** 是一种基于荧光的SNP基因分型技术。
    
    ### 输入格式
    - SNP位点使用方括号标记: `[A/G]`, `[C/T]` 等
    - 上游序列建议 ≥50bp
    - 下游序列建议 ≥200bp
    
    ### 输出说明
    - **FAM引物**: 检测第一个等位基因
    - **HEX引物**: 检测第二个等位基因  
    - **通用引物**: 反向引物，两个等位基因共用
    
    ### 质量评分
    | 等级 | 分数范围 | 建议 |
    |------|----------|------|
    | 优秀 ★★★★★ | ≥85 | 可直接使用 |
    | 良好 ★★★★☆ | 75-84 | 推荐使用 |
    | 合格 ★★★☆☆ | 65-74 | 可以使用 |
    | 可用 ★★☆☆☆ | 55-64 | 需要优化 |
    | 勉强 ★☆☆☆☆ | 45-54 | 建议调整 |
    | 需优化 ☆☆☆☆☆ | <45 | 不建议使用 |
    
    ---
    
    ## 常规PCR引物设计
    
    用于设计普通PCR扩增引物对。
    
    ### 设计原则
    - 引物长度: 18-25 bp
    - Tm值: 55-68°C
    - GC含量: 40-60%
    - 产物大小: 可自定义
    
    ### 注意事项
    - 避免引物3'端自身互补
    - 避免引物对之间形成二聚体
    - 两条引物Tm差异应 <2°C
    
    ---
    
    ## 常见问题
    
    **Q: 为什么设计不出引物？**
    A: 可能原因：
    1. 序列太短
    2. GC含量过高或过低
    3. 序列中有过多重复
    
    **Q: 如何提高引物质量？**
    A: 
    1. 提供更长的侧翼序列
    2. 选择GC含量适中的区域
    3. 调整参数设置
    """)


# ==================== 主程序 ====================

def main():
    # 侧边栏导航
    st.sidebar.markdown("## 🧬 引物设计工具")
    st.sidebar.markdown("**v5.0 Web版**")
    st.sidebar.markdown("---")
    
    page = st.sidebar.radio(
        "选择功能",
        ["🏠 首页", "🔬 KASP引物设计", "🧪 常规PCR引物设计", "🔍 引物分析", "📖 帮助文档"]
    )
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    <small>
    
    **关于本工具**
    
    本工具用于设计KASP基因分型引物和常规PCR引物。
    
    支持功能:
    - KASP引物多方案设计
    - 常规PCR引物对设计
    - 引物质量综合评估
    - CSV报告导出
    
    </small>
    """, unsafe_allow_html=True)
    
    # 页面路由
    if page == "🏠 首页":
        st.markdown('<p class="main-header">🧬 引物设计工具</p>', unsafe_allow_html=True)
        st.markdown('<p class="sub-header">KASP & 常规PCR 引物设计平台 v5.0</p>', unsafe_allow_html=True)
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            ### 🔬 KASP引物设计
            
            针对SNP位点设计KASP基因分型引物
            
            - ✅ 自动添加FAM/HEX荧光尾巴
            - ✅ 智能错配位点设计
            - ✅ 多方案评分排序
            - ✅ 一键导出CSV
            
            """)
            if st.button("开始KASP设计 →", key="goto_kasp"):
                st.session_state['page'] = "🔬 KASP引物设计"
                st.rerun()
        
        with col2:
            st.markdown("""
            ### 🧪 常规PCR引物设计
            
            设计普通PCR扩增引物对
            
            - ✅ 自定义产物大小
            - ✅ 指定目标区域
            - ✅ 引物对质量评估
            - ✅ 二聚体风险检测
            
            """)
            if st.button("开始PCR设计 →", key="goto_pcr"):
                st.session_state['page'] = "🧪 常规PCR引物设计"
                st.rerun()
        
        st.markdown("---")
        
        st.markdown("""
        ### 🔍 快速引物分析
        
        输入任意引物序列，快速评估其质量指标
        """)
        
        quick_primer = st.text_input("输入引物序列进行快速分析", placeholder="例如: ATGCGATCGATCGATCGATCG")
        
        if quick_primer:
            primer = re.sub(r'[^ATGC]', '', quick_primer.upper())
            if len(primer) >= 10:
                result = evaluate_primer_quality(primer)
                grade, stars, _ = get_quality_grade(result['score'])
                
                c1, c2, c3, c4, c5 = st.columns(5)
                c1.metric("长度", f"{result['length']}bp")
                c2.metric("Tm", f"{result['tm']}°C")
                c3.metric("GC%", f"{result['gc_content']:.1f}%")
                c4.metric("评分", f"{result['score']:.0f}")
                c5.metric("等级", f"{grade}")
    
    elif page == "🔬 KASP引物设计":
        show_kasp_design()
    
    elif page == "🧪 常规PCR引物设计":
        show_regular_pcr_design()
    
    elif page == "🔍 引物分析":
        show_primer_analysis()
    
    elif page == "📖 帮助文档":
        show_help()


if __name__ == "__main__":
    main()
