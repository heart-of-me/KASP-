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
    MIN_GC: float = 30.0  # 小麦大忌：避免<30%
    MAX_GC: float = 65.0  # 小麦大忌：避免>65%
    OPTIMAL_GC_MIN: float = 40.0
    OPTIMAL_GC_MAX: float = 55.0
    MAX_TM_DIFF: float = 2.0
    # KASP产物大小：小麦建议50-100bp
    REV_MIN_DISTANCE: int = 30   # 最近距离，产物约50bp
    REV_MAX_DISTANCE: int = 80   # 最远距离，产物约100bp
    PRODUCT_MIN: int = 50        # 小麦大忌#3：扩增子过长
    PRODUCT_MAX: int = 120       # KASP建议50-100bp
    MISMATCH_POSITIONS: List[int] = None
    
    # 🌾 小麦特异性参数
    WHEAT_MODE: bool = False
    WHEAT_CHECK_FLANKING_SNP: bool = True   # 大忌#2：检查侧翼干扰SNP
    WHEAT_CHECK_REPEAT: bool = True         # 大忌#5：检查重复序列
    WHEAT_STRICT_GC: bool = True            # 大忌#4：GC含量严格模式
    
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
    # 小麦特异性参数
    WHEAT_MODE: bool = False
    WHEAT_AVOID_5PRIME_PERCENT: float = 40.0  # 避开5'端的百分比
    WHEAT_PREFER_3PRIME: bool = True  # 优先3'端区域
    WHEAT_MIN_UNIQUE_BASES: int = 3  # 3'端最少独特碱基数


# ==================== 核心计算函数 ====================

def calc_gc_content(seq: str) -> float:
    """计算GC含量百分比"""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def calc_tm_nearest_neighbor(seq: str, na_conc: float = 50.0, primer_conc: float = 250.0) -> float:
    """
    使用最近邻法(Nearest-Neighbor)计算Tm值
    这是目前最准确的Tm计算方法，与实际实验值偏差通常<2°C
    
    参数:
        seq: 引物序列
        na_conc: Na+浓度 (mM)，默认50mM
        primer_conc: 引物浓度 (nM)，默认250nM
    
    参考: SantaLucia J Jr. (1998) PNAS 95:1460-1465
    """
    seq = seq.upper()
    length = len(seq)
    
    if length < 8:
        # 序列过短，使用简单公式
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        return round(2 * at_count + 4 * gc_count, 1)
    
    # 最近邻热力学参数 (ΔH kcal/mol, ΔS cal/mol·K)
    # SantaLucia 1998 统一参数
    nn_params = {
        'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
        'AT': (-7.2, -20.4),
        'TA': (-7.2, -21.3),
        'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
        'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
        'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
        'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
        'CG': (-10.6, -27.2),
        'GC': (-9.8, -24.4),
        'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
    }
    
    # 起始参数
    init_params = {
        'G': (0.1, -2.8), 'C': (0.1, -2.8),
        'A': (2.3, 4.1), 'T': (2.3, 4.1)
    }
    
    # 计算ΔH和ΔS
    dH = 0.0  # kcal/mol
    dS = 0.0  # cal/mol·K
    
    # 起始贡献
    if seq[0] in init_params:
        dH += init_params[seq[0]][0]
        dS += init_params[seq[0]][1]
    if seq[-1] in init_params:
        dH += init_params[seq[-1]][0]
        dS += init_params[seq[-1]][1]
    
    # 最近邻贡献
    for i in range(length - 1):
        dinuc = seq[i:i+2]
        if dinuc in nn_params:
            dH += nn_params[dinuc][0]
            dS += nn_params[dinuc][1]
    
    # 盐浓度校正 (von Ahsen 2001)
    # ΔS_corrected = ΔS + 0.368 * N * ln([Na+])
    import math
    na_molar = na_conc / 1000.0  # 转换为M
    dS_corrected = dS + 0.368 * (length - 1) * math.log(na_molar)
    
    # 计算Tm
    # Tm = ΔH / (ΔS + R * ln(Ct/4)) - 273.15
    # R = 1.987 cal/mol·K
    R = 1.987
    ct = primer_conc * 1e-9  # 转换为M
    
    tm = (dH * 1000) / (dS_corrected + R * math.log(ct / 4)) - 273.15
    
    return round(tm, 1)


def calc_tm_simple(seq: str) -> float:
    """
    使用改进的Wallace公式计算Tm值（用于快速筛选）
    适用于14-30bp的引物
    """
    seq = seq.upper()
    length = len(seq)
    gc_count = seq.count('G') + seq.count('C')
    
    if length < 14:
        # Wallace公式
        tm = 2 * (length - gc_count) + 4 * gc_count
    else:
        # 改进的公式，考虑长度影响
        gc_percent = (gc_count / length) * 100
        # Primer3使用的公式变体
        tm = 64.9 + 41 * (gc_count - 16.4) / length
    
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

# ========== 🌾 小麦KASP五大忌检测函数 ==========

def check_flanking_snp_risk(upstream: str, downstream: str, primer_region_len: int = 25) -> Tuple[bool, List[str]]:
    """
    大忌#2：检查SNP侧翼干扰
    检测引物结合区域内是否有潜在的SNP/InDel风险
    
    返回: (is_risky, issues_list)
    """
    issues = []
    is_risky = False
    
    # 检查上游区域（allele-specific引物区域）
    upstream_region = upstream[-primer_region_len:] if len(upstream) >= primer_region_len else upstream
    
    # 检测低复杂度区域（可能是潜在变异热点）
    for base in 'ATGC':
        if base * 4 in upstream_region:
            issues.append(f"Allele引物区域含{base}x4+重复，可能存在InDel变异")
            is_risky = True
    
    # 检测二核苷酸重复（常见SNP位点）
    for i in range(len(upstream_region) - 5):
        dinuc = upstream_region[i:i+2]
        if upstream_region.count(dinuc) >= 4:
            issues.append(f"Allele引物区域含{dinuc}重复≥×4，SNP风险高")
            is_risky = True
            break
    
    # 检查下游区域（Common引物区域）
    # Common引物通常在30-100bp处
    if len(downstream) >= 30:
        common_region = downstream[30:80] if len(downstream) >= 80 else downstream[30:]
        for base in 'ATGC':
            if base * 4 in common_region:
                issues.append(f"Common引物区域含{base}x4+重复")
    
    return is_risky, issues


def check_wheat_repeat_sequences(seq: str) -> Tuple[bool, List[str]]:
    """
    大忌#5：检查小麦重复序列
    检测常见的小麦转座子/重复序列特征
    """
    seq = seq.upper()
    issues = []
    has_repeat = False
    
    # 小麦常见转座子特征序列（简化版）
    repeat_motifs = [
        ('CACTA', 'CACTA转座子'),
        ('CACTG', 'CACTA类转座子'),
        ('TGCA', 'Mariner类转座子'),
        ('GGCCGG', '稀有切割位点/转座子'),
        ('ATATATAT', 'AT-rich重复'),
        ('GCGCGCGC', 'GC-rich重复'),
        ('CCCGGG', 'SmaI位点/重复区'),
    ]
    
    for motif, desc in repeat_motifs:
        if motif in seq:
            issues.append(f"含{desc}特征")
            has_repeat = True
    
    # 检测简单重复序列(SSR)
    # 二核苷酸重复
    for dinuc in ['AT', 'TA', 'GC', 'CG', 'AG', 'GA', 'CT', 'TC', 'AC', 'CA', 'GT', 'TG']:
        if (dinuc * 5) in seq:  # 5次以上重复
            issues.append(f"({dinuc})n SSR重复")
            has_repeat = True
            break
    
    # 三核苷酸重复
    common_trinuc = ['AAT', 'ATT', 'GGC', 'GCC', 'GAA', 'TTC']
    for trinuc in common_trinuc:
        if (trinuc * 4) in seq:
            issues.append(f"({trinuc})n SSR重复")
            has_repeat = True
            break
    
    return has_repeat, issues


def check_gc_extreme(seq: str, strict: bool = True) -> Tuple[bool, float, str]:
    """
    大忌#4：检查GC含量极端
    小麦KASP严格模式：30-65%
    
    返回: (is_extreme, gc_content, message)
    """
    gc = calc_gc_content(seq)
    
    if strict:
        # 小麦严格模式
        if gc > 65:
            return True, gc, f"GC过高({gc:.1f}%)>小麦大忌65%，易形成二级结构"
        elif gc < 30:
            return True, gc, f"GC过低({gc:.1f}%)<小麦大忌30%，Tm不足"
        elif gc > 60:
            return False, gc, f"GC偏高({gc:.1f}%)，注意二级结构"
        elif gc < 35:
            return False, gc, f"GC偏低({gc:.1f}%)，注意Tm值"
    else:
        if gc > 70 or gc < 25:
            return True, gc, f"GC极端({gc:.1f}%)"
    
    return False, gc, "GC含量正常"


def check_amplicon_length_kasp(product_size: int) -> Tuple[str, int]:
    """
    大忌#3：检查KASP扩增子长度
    KASP建议50-100bp
    
    返回: (评价, 评分加成)
    """
    if 50 <= product_size <= 80:
        return "★ 最佳(50-80bp)", 15
    elif 80 < product_size <= 100:
        return "✓ 良好(80-100bp)", 10
    elif 100 < product_size <= 120:
        return "△ 可接受(100-120bp)", 0
    elif 120 < product_size <= 150:
        return "⚠ 偏长(120-150bp)", -10
    else:
        return "✗ 过长(>150bp)，小麦大忌", -25


def evaluate_kasp_wheat_specificity(upstream: str, downstream: str, 
                                     fwd_seq: str, rev_seq: str,
                                     config: KASPConfig) -> Tuple[float, List[str], Dict]:
    """
    综合评估小麦KASP引物的特异性
    
    返回: (评分加成, 问题列表, 详细信息)
    """
    score_bonus = 0
    issues = []
    details = {}
    
    # 1. 检查侧翼干扰（大忌#2）
    if config.WHEAT_CHECK_FLANKING_SNP:
        flanking_risky, flanking_issues = check_flanking_snp_risk(upstream, downstream)
        if flanking_risky:
            score_bonus -= 15
            issues.extend(flanking_issues)
        details['flanking_risk'] = flanking_risky
    
    # 2. 检查重复序列（大忌#5）
    if config.WHEAT_CHECK_REPEAT:
        # 检查整个区域
        full_seq = upstream[-30:] + downstream[:100]
        has_repeat, repeat_issues = check_wheat_repeat_sequences(full_seq)
        if has_repeat:
            score_bonus -= 20
            issues.extend(["🌾 " + i for i in repeat_issues])
        details['has_repeat'] = has_repeat
        
        # 单独检查引物
        fwd_repeat, fwd_rep_issues = check_wheat_repeat_sequences(fwd_seq)
        rev_repeat, rev_rep_issues = check_wheat_repeat_sequences(rev_seq)
        if fwd_repeat or rev_repeat:
            score_bonus -= 10
            issues.append("引物序列含重复元件")
    
    # 3. GC含量检查（大忌#4）
    if config.WHEAT_STRICT_GC:
        fwd_gc_extreme, fwd_gc, fwd_gc_msg = check_gc_extreme(fwd_seq, strict=True)
        rev_gc_extreme, rev_gc, rev_gc_msg = check_gc_extreme(rev_seq, strict=True)
        
        if fwd_gc_extreme:
            score_bonus -= 20
            issues.append(f"Allele引物{fwd_gc_msg}")
        if rev_gc_extreme:
            score_bonus -= 15
            issues.append(f"Common引物{rev_gc_msg}")
        
        details['fwd_gc'] = fwd_gc
        details['rev_gc'] = rev_gc
    
    # 4. 3'端特异性检查（小麦同源基因关键）
    fwd_end_3 = fwd_seq[-3:].upper()
    # 小麦中常见的4fdd保守3'端
    conserved_ends = ['GGC', 'GCC', 'CGG', 'CCG', 'GCG', 'CGC']
    if fwd_end_3 in conserved_ends:
        score_bonus -= 5
        issues.append(f"Allele引灩3'端{fwd_end_3}在小麦中较保守")
    
    # 5. 序列复杂度检查
    complexity = analyze_sequence_complexity(fwd_seq)
    if complexity['complexity_score'] < 50:
        score_bonus -= 10
        issues.append("引物序列复杂度低，可能匹配多个位点")
    elif complexity['complexity_score'] > 75:
        score_bonus += 5
    
    details['complexity_score'] = complexity['complexity_score']
    
    return score_bonus, issues, details


def generate_wheat_warnings(config: KASPConfig) -> List[str]:
    """生成小麦模式警告信息"""
    warnings = []
    warnings.append("🌾 小麦模式已启用 - 请注意以下五大忌：")
    warnings.append("1️⃣ 大忌#1 同源基因干扰：请将引物BLAST到A/B/D三个基因组验证特异性")
    warnings.append("2️⃣ 大忌#2 侧翼干扰：确认Allele引物区域无其他SNP/InDel")
    warnings.append("3️⃣ 大忌#3 扩增子过长：已优化为50-100bp")
    warnings.append("4️⃣ 大忌#4 GC极端：已检测30-65%范围")
    warnings.append("5️⃣ 大忌#5 重复序列：已检测转座子/SSR特征")
    return warnings


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
    """
    设计多套KASP引物方案 - 优化版
    确保不产生重复引物，质量不达标时返回空列表
    """
    if config is None:
        config = KASPConfig()
    
    # 检查序列长度是否足够
    if len(upstream) < config.MIN_PRIMER_LEN:
        return []  # 上游序列太短
    
    if len(downstream) < config.REV_MIN_DISTANCE + config.MIN_PRIMER_LEN:
        return []  # 下游序列太短
    
    all_schemes = []
    
    # 生成不同长度的正向引物
    for primer_len in range(config.MIN_PRIMER_LEN, min(config.MAX_PRIMER_LEN + 1, len(upstream) + 1)):
        core_seq = upstream[-(primer_len - 1):]
        
        if len(core_seq) < config.MIN_PRIMER_LEN - 1:
            continue
        
        # 小麦模式：预筛选GC含量
        if config.WHEAT_MODE and config.WHEAT_STRICT_GC:
            core_gc = calc_gc_content(core_seq)
            if core_gc > 65 or core_gc < 30:
                continue  # 跳过GC极端的候选
        
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
            
            # 质量过低的直接跳过
            if eval1['score'] < 40 or eval2['score'] < 40:
                continue
            
            tm_diff = abs(eval1['tm'] - eval2['tm'])
            
            # Tm差异过大直接跳过
            if tm_diff > config.MAX_TM_DIFF + 1:
                continue
            
            # 搜索反向引物
            max_rev_dist = min(config.REV_MAX_DISTANCE + 1, len(downstream) - config.MIN_PRIMER_LEN + 1)
            if max_rev_dist <= config.REV_MIN_DISTANCE:
                continue
                
            for rev_dist in range(config.REV_MIN_DISTANCE, max_rev_dist):
                for rev_len in range(config.MIN_PRIMER_LEN, min(config.MAX_PRIMER_LEN + 1, len(downstream) - rev_dist + 1)):
                    rev_start = rev_dist
                    rev_end = rev_dist + rev_len
                    
                    if rev_end > len(downstream):
                        continue
                    
                    rev_seq = reverse_complement(downstream[rev_start:rev_end])
                    
                    # 小麦模式：预筛选反向引物GC
                    if config.WHEAT_MODE and config.WHEAT_STRICT_GC:
                        rev_gc = calc_gc_content(rev_seq)
                        if rev_gc > 65 or rev_gc < 30:
                            continue
                    
                    eval_rev = evaluate_primer_quality(rev_seq, config)
                    
                    # 质量过低跳过
                    if eval_rev['score'] < 40:
                        continue
                    
                    # 检查引物二聚体
                    has_dimer = (check_primer_dimer(fwd_allele1, rev_seq) or 
                                check_primer_dimer(fwd_allele2, rev_seq))
                    
                    # 计算产物大小
                    product_size = len(upstream) + 1 + rev_dist + rev_len
                    
                    # === 基础评分 ===
                    avg_fwd_score = (eval1['score'] + eval2['score']) / 2
                    total_score = (avg_fwd_score * 0.4 + eval_rev['score'] * 0.3)
                    
                    # Tm匹配评分
                    if tm_diff <= 0.5:
                        total_score += 15
                    elif tm_diff <= 1.0:
                        total_score += 10
                    elif tm_diff <= 2.0:
                        total_score += 5
                    else:
                        total_score -= 10
                    
                    # 二聚体惩罚
                    if has_dimer:
                        total_score -= 15
                    
                    # === 小麦特异性评分（五大忌）===
                    wheat_issues = []
                    wheat_details = {}
                    
                    if config.WHEAT_MODE:
                        wheat_bonus, wheat_issues, wheat_details = evaluate_kasp_wheat_specificity(
                            upstream, downstream, fwd_allele1, rev_seq, config
                        )
                        total_score += wheat_bonus
                        
                        # 大忌#3：扩增子长度评分
                        amplicon_status, amplicon_bonus = check_amplicon_length_kasp(product_size)
                        total_score += amplicon_bonus
                        wheat_details['amplicon_status'] = amplicon_status
                    else:
                        # 非小麦模式的产物大小评分
                        if config.PRODUCT_MIN <= product_size <= config.PRODUCT_MAX:
                            total_score += 5
                    
                    total_score = max(0, min(100, total_score))
                    
                    # 判断是否可用（小麦模式更严格）
                    is_usable = True
                    if config.WHEAT_MODE:
                        is_usable = (
                            total_score >= 50 and
                            product_size <= 120 and
                            30 <= eval1['gc_content'] <= 65 and
                            30 <= eval_rev['gc_content'] <= 65 and
                            not has_dimer and
                            tm_diff <= config.MAX_TM_DIFF
                        )
                    else:
                        is_usable = (
                            total_score >= 45 and
                            not has_dimer and
                            tm_diff <= config.MAX_TM_DIFF
                        )
                    
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
                        'total_score': total_score,
                        'is_usable': is_usable,
                        'wheat_mode': config.WHEAT_MODE,
                        'wheat_issues': wheat_issues,
                        'wheat_details': wheat_details
                    }
                    all_schemes.append(scheme)
    
    # 如果没有找到任何方案，返回空列表
    if not all_schemes:
        return []
    
    # 按评分排序
    if config.WHEAT_MODE:
        all_schemes.sort(key=lambda x: (x.get('is_usable', False), x['total_score']), reverse=True)
    else:
        all_schemes.sort(key=lambda x: x['total_score'], reverse=True)
    
    # 严格去重：确保每个方案的引物组合都是唯一的
    unique_schemes = []
    seen_cores = set()  # 核心序列组合
    seen_full = set()   # 完整引物组合
    
    for scheme in all_schemes:
        # 多重去重检查
        core_key = (scheme['fwd_allele1_core'], scheme['reverse'])
        full_key = (scheme['fwd_allele1_full'], scheme['fwd_allele2_full'], scheme['reverse'])
        
        if core_key in seen_cores or full_key in seen_full:
            continue
        
        seen_cores.add(core_key)
        seen_full.add(full_key)
        unique_schemes.append(scheme)
        
        if len(unique_schemes) >= num_schemes:
            break
    
    return unique_schemes
        if key not in seen:
            seen.add(key)
            unique_schemes.append(scheme)
            if len(unique_schemes) >= num_schemes:
                break
    
    return unique_schemes


# ==================== 常规PCR引物设计 ====================

def analyze_sequence_complexity(seq: str) -> Dict:
    """
    分析序列复杂度，用于评估小麦同源基因特异性
    返回各种复杂度指标
    """
    seq = seq.upper()
    length = len(seq)
    
    # 计算不同碱基的分布
    base_counts = {base: seq.count(base) for base in 'ATGC'}
    
    # 计算二核苷酸多样性
    dinucs = [seq[i:i+2] for i in range(len(seq)-1)]
    unique_dinucs = len(set(dinucs))
    dinuc_diversity = unique_dinucs / 16  # 最多16种二核苷酸
    
    # 计算三核苷酸多样性（密码子相关）
    trinucs = [seq[i:i+3] for i in range(len(seq)-2)]
    unique_trinucs = len(set(trinucs))
    trinuc_diversity = unique_trinucs / min(64, len(trinucs)) if trinucs else 0
    
    # 检测低复杂度区域
    low_complexity = False
    for base in 'ATGC':
        if base * 4 in seq:  # 4个连续相同碱基
            low_complexity = True
            break
    
    # 整体复杂度评分 (0-100)
    complexity_score = (dinuc_diversity * 50 + trinuc_diversity * 50)
    
    return {
        'length': length,
        'base_counts': base_counts,
        'dinuc_diversity': dinuc_diversity,
        'trinuc_diversity': trinuc_diversity,
        'low_complexity': low_complexity,
        'complexity_score': complexity_score
    }


def check_wheat_specificity(seq: str, position_in_gene: float) -> Tuple[float, List[str]]:
    """
    评估引物在小麦中的特异性潜力
    
    参数:
        seq: 引物序列
        position_in_gene: 引物在基因中的相对位置 (0-1, 0=5'端, 1=3'端)
    
    返回:
        (特异性评分加成, 问题列表)
    """
    seq = seq.upper()
    issues = []
    score_bonus = 0
    
    # 1. 位置评估 - 5'端保守区域风险
    if position_in_gene < 0.2:
        issues.append("位于5'端高度保守区(可能扩增A/B/D三拷贝)")
        score_bonus -= 25
    elif position_in_gene < 0.4:
        issues.append("位于5'端较保守区")
        score_bonus -= 10
    elif position_in_gene > 0.7:
        score_bonus += 10  # 3'端/UTR区域，通常变异更多
    
    # 2. 3'端碱基多样性检查（关键！）
    end_6 = seq[-6:]
    unique_end_bases = len(set(end_6))
    if unique_end_bases >= 4:
        score_bonus += 8  # 3'端碱基多样，有利于特异性
    elif unique_end_bases <= 2:
        issues.append("3'端碱基单一，特异性可能不足")
        score_bonus -= 8
    
    # 3. 检查是否含有小麦中常见的保守基序
    # 常见的起始密码子附近保守序列
    conserved_motifs = [
        'ATGGC', 'ATGGA', 'ATGGG',  # 起始密码子附近
        'GCGGC', 'GCCGC', 'CGGCG',  # 高GC保守区
        'AAGAAG', 'GAAGAA',  # 常见重复
    ]
    for motif in conserved_motifs:
        if motif in seq:
            issues.append(f"含保守基序{motif}")
            score_bonus -= 5
            break
    
    # 4. 检查3'端最后3个碱基的独特性
    end_3 = seq[-3:]
    # 避免常见的保守3'端
    common_ends = ['GGC', 'GCC', 'CGG', 'CCG', 'AAA', 'TTT']
    if end_3 in common_ends:
        issues.append(f"3'端{end_3}较常见")
        score_bonus -= 5
    
    # 5. 序列复杂度
    complexity = analyze_sequence_complexity(seq)
    if complexity['complexity_score'] < 60:
        issues.append("序列复杂度较低")
        score_bonus -= 5
    elif complexity['complexity_score'] > 80:
        score_bonus += 5
    
    return score_bonus, issues


def check_repeat_region(seq: str, max_repeat: int = 4) -> bool:
    """检查是否有连续重复碱基"""
    seq = seq.upper()
    for base in ['A', 'T', 'G', 'C']:
        if base * max_repeat in seq:
            return True
    # 检查二核苷酸重复
    for i in range(len(seq) - 5):
        dinuc = seq[i:i+2]
        if dinuc * 3 in seq:
            return True
    return False


def check_gc_clamp(seq: str) -> Tuple[bool, str]:
    """检查GC夹（3'端最后5个碱基中的GC数量）"""
    end_5 = seq[-5:].upper()
    gc_count = end_5.count('G') + end_5.count('C')
    last_2 = seq[-2:].upper()
    last_2_gc = last_2.count('G') + last_2.count('C')
    
    # 理想情况：3'端有1-2个GC
    if 1 <= last_2_gc <= 2 and 2 <= gc_count <= 3:
        return True, "GC夹良好"
    elif gc_count > 3:
        return False, "3'端GC过多，可能非特异性扩增"
    elif gc_count == 0:
        return False, "3'端无GC，结合力弱"
    else:
        return True, "GC夹可接受"


def evaluate_primer_quality_strict(seq: str, config) -> Dict:
    """严格评估引物质量（用于常规PCR）"""
    seq = seq.upper()
    
    gc_content = calc_gc_content(seq)
    tm = calc_tm_nearest_neighbor(seq)
    has_hairpin = check_hairpin(seq)
    has_self_dimer = check_self_dimer(seq)
    three_prime_ok, three_prime_msg = check_gc_clamp(seq)
    has_repeat = check_repeat_region(seq)
    
    result = {
        'sequence': seq,
        'length': len(seq),
        'gc_content': gc_content,
        'tm': tm,
        'has_hairpin': has_hairpin,
        'has_self_dimer': has_self_dimer,
        'has_repeat': has_repeat,
        'three_prime_ok': three_prime_ok,
        'three_prime_msg': three_prime_msg,
        'issues': [],
        'score': 100
    }
    
    # === Tm评分（权重高）===
    if tm < config.MIN_TM:
        result['issues'].append(f"Tm过低({tm}°C < {config.MIN_TM}°C)")
        result['score'] -= 20
    elif tm > config.MAX_TM:
        result['issues'].append(f"Tm过高({tm}°C > {config.MAX_TM}°C)")
        result['score'] -= 15
    elif abs(tm - config.OPTIMAL_TM) <= 2:
        result['score'] += 5  # 接近最优Tm加分
    elif abs(tm - config.OPTIMAL_TM) > 5:
        result['issues'].append(f"Tm偏离最优值({tm}°C)")
        result['score'] -= 5
    
    # === GC含量评分 ===
    if gc_content < config.MIN_GC:
        result['issues'].append(f"GC含量过低({gc_content:.1f}% < {config.MIN_GC}%)")
        result['score'] -= 20
    elif gc_content > config.MAX_GC:
        result['issues'].append(f"GC含量过高({gc_content:.1f}% > {config.MAX_GC}%)")
        result['score'] -= 20
    elif config.OPTIMAL_GC_MIN <= gc_content <= config.OPTIMAL_GC_MAX:
        result['score'] += 5  # 最优GC范围加分
    
    # === 二级结构评分 ===
    if has_hairpin:
        result['issues'].append("可能形成发夹结构")
        result['score'] -= 15
    
    if has_self_dimer:
        result['issues'].append("自身二聚体风险")
        result['score'] -= 15
    
    # === 重复序列检测 ===
    if has_repeat:
        result['issues'].append("含有连续重复碱基")
        result['score'] -= 10
    
    # === 3'端稳定性 ===
    if not three_prime_ok:
        result['issues'].append(three_prime_msg)
        result['score'] -= 15
    
    # === 长度评分 ===
    if len(seq) < config.MIN_PRIMER_LEN:
        result['issues'].append(f"长度过短({len(seq)}bp)")
        result['score'] -= 15
    elif len(seq) > config.MAX_PRIMER_LEN:
        result['issues'].append(f"长度过长({len(seq)}bp)")
        result['score'] -= 10
    elif config.MIN_PRIMER_LEN + 2 <= len(seq) <= config.OPTIMAL_PRIMER_LEN + 2:
        result['score'] += 3  # 最优长度加分
    
    # === 5'端检查 ===
    if seq[0] in ['G', 'C']:
        result['score'] += 2  # 5'端GC加分
    
    result['score'] = max(0, min(100, result['score']))
    
    return result


def design_regular_primers(sequence: str, config: RegularPCRConfig = None, 
                          num_pairs: int = 5, target_start: int = None, 
                          target_end: int = None) -> List[Dict]:
    """设计常规PCR引物对 - 优化版（支持小麦模式）"""
    if config is None:
        config = RegularPCRConfig()
    
    sequence = re.sub(r'[^ATGC]', '', sequence.upper())
    seq_len = len(sequence)
    
    # === 小麦模式：避开5'端保守区 ===
    wheat_avoid_region = 0
    if config.WHEAT_MODE:
        # 计算需要避开的5'端区域长度
        wheat_avoid_region = int(seq_len * config.WHEAT_AVOID_5PRIME_PERCENT / 100)
        wheat_avoid_region = max(wheat_avoid_region, 50)  # 至少避开50bp
    
    # 确定目标区域
    if target_start is None and target_end is None:
        # 未指定目标区域时
        margin = max(config.MIN_PRIMER_LEN + 10, 30)
        
        if config.WHEAT_MODE:
            # 小麦模式：从避开区域之后开始
            target_start = max(margin, wheat_avoid_region)
            target_end = seq_len - margin
        else:
            target_start = margin
            target_end = seq_len - margin
    else:
        if target_start is None:
            target_start = max(30, target_end - config.PRODUCT_MAX)
        if target_end is None:
            target_end = min(seq_len - 30, target_start + config.PRODUCT_MAX)
        
        # 小麦模式下，警告用户如果目标区域在5'端
        if config.WHEAT_MODE and target_start < wheat_avoid_region:
            # 调整目标起始位置
            target_start = max(target_start, wheat_avoid_region)
    
    all_pairs = []
    
    # 计算有效的搜索范围
    if config.WHEAT_MODE:
        # 小麦模式：正向引物从避开区域之后开始搜索
        fwd_search_start = max(wheat_avoid_region, target_start - config.PRODUCT_MAX // 2)
        # 优先在3'端区域设计
        if config.WHEAT_PREFER_3PRIME:
            preferred_start = int(seq_len * 0.4)  # 从40%位置开始
            fwd_search_start = max(fwd_search_start, preferred_start)
    else:
        fwd_search_start = max(0, target_start - config.PRODUCT_MAX)
    
    fwd_search_end = max(fwd_search_start + 20, target_start - config.MIN_PRIMER_LEN)
    fwd_search_end = min(fwd_search_end, seq_len - config.PRODUCT_MIN - config.MIN_PRIMER_LEN)
    
    # 反向引物搜索区域
    rev_search_start = min(target_end + config.MIN_PRIMER_LEN, seq_len - config.MIN_PRIMER_LEN)
    rev_search_start = max(rev_search_start, config.PRODUCT_MIN + config.MIN_PRIMER_LEN)
    rev_search_end = min(target_end + config.PRODUCT_MAX, seq_len)
    
    # 如果搜索范围无效，扩展搜索
    if fwd_search_end <= fwd_search_start:
        if config.WHEAT_MODE:
            fwd_search_start = wheat_avoid_region
            fwd_search_end = max(fwd_search_start + 50, seq_len * 2 // 3)
        else:
            fwd_search_start = 0
            fwd_search_end = max(50, seq_len // 3)
    
    if rev_search_end <= rev_search_start:
        rev_search_start = max(seq_len * 2 // 3, config.PRODUCT_MIN)
        rev_search_end = seq_len
    
    # 使用更智能的步进策略
    step = max(3, (fwd_search_end - fwd_search_start) // 30)
    
    for fwd_start in range(fwd_search_start, fwd_search_end, step):
        for fwd_len in range(config.MIN_PRIMER_LEN, config.MAX_PRIMER_LEN + 1):
            fwd_end = fwd_start + fwd_len
            if fwd_end > seq_len:
                continue
            
            fwd_seq = sequence[fwd_start:fwd_end]
            
            # 快速预筛选
            gc = calc_gc_content(fwd_seq)
            if gc < config.MIN_GC - 5 or gc > config.MAX_GC + 5:
                continue
            
            fwd_eval = evaluate_primer_quality_strict(fwd_seq, config)
            
            # 严格的质量门槛
            if fwd_eval['score'] < 50:
                continue
            
            # 搜索匹配的反向引物
            rev_step = max(3, (rev_search_end - rev_search_start) // 30)
            
            for rev_end in range(rev_search_start, rev_search_end, rev_step):
                for rev_len in range(config.MIN_PRIMER_LEN, config.MAX_PRIMER_LEN + 1):
                    rev_start = rev_end - rev_len
                    if rev_start < 0 or rev_start <= fwd_end:
                        continue
                    
                    product_size = rev_end - fwd_start
                    if not (config.PRODUCT_MIN <= product_size <= config.PRODUCT_MAX):
                        continue
                    
                    rev_seq = reverse_complement(sequence[rev_start:rev_end])
                    
                    # 快速预筛选
                    gc = calc_gc_content(rev_seq)
                    if gc < config.MIN_GC - 5 or gc > config.MAX_GC + 5:
                        continue
                    
                    rev_eval = evaluate_primer_quality_strict(rev_seq, config)
                    
                    if rev_eval['score'] < 50:
                        continue
                    
                    tm_diff = abs(fwd_eval['tm'] - rev_eval['tm'])
                    
                    # Tm差异过大直接跳过
                    if tm_diff > config.MAX_TM_DIFF + 1:
                        continue
                    
                    has_dimer = check_primer_dimer(fwd_seq, rev_seq)
                    
                    # === 综合评分计算 ===
                    # 基础分：两个引物评分的平均值
                    base_score = (fwd_eval['score'] + rev_eval['score']) / 2
                    
                    # Tm匹配奖惩（重要指标）
                    if tm_diff <= 0.5:
                        tm_bonus = 12
                    elif tm_diff <= 1.0:
                        tm_bonus = 8
                    elif tm_diff <= 1.5:
                        tm_bonus = 4
                    elif tm_diff <= 2.0:
                        tm_bonus = 0
                    else:
                        tm_bonus = -15
                    
                    # 二聚体惩罚
                    dimer_penalty = -20 if has_dimer else 0
                    
                    # 产物大小奖励（最优范围200-400bp）
                    if 200 <= product_size <= 400:
                        size_bonus = 5
                    elif 150 <= product_size <= 500:
                        size_bonus = 2
                    else:
                        size_bonus = 0
                    
                    # Tm绝对值检查（确保都在合理范围）
                    avg_tm = (fwd_eval['tm'] + rev_eval['tm']) / 2
                    if 58 <= avg_tm <= 62:
                        tm_range_bonus = 5
                    elif 55 <= avg_tm <= 65:
                        tm_range_bonus = 2
                    else:
                        tm_range_bonus = -5
                    
                    # === 小麦特异性评分 ===
                    wheat_bonus = 0
                    wheat_issues = []
                    if config.WHEAT_MODE:
                        # 计算引物在序列中的相对位置
                        fwd_position = fwd_start / seq_len
                        rev_position = rev_end / seq_len
                        
                        # 评估正向引物特异性
                        fwd_wheat_bonus, fwd_wheat_issues = check_wheat_specificity(fwd_seq, fwd_position)
                        # 评估反向引物特异性
                        rev_wheat_bonus, rev_wheat_issues = check_wheat_specificity(rev_seq, rev_position)
                        
                        wheat_bonus = (fwd_wheat_bonus + rev_wheat_bonus) / 2
                        wheat_issues = fwd_wheat_issues + rev_wheat_issues
                        
                        # 额外奖励：两个引物都在3'半区
                        if fwd_position > 0.5 and rev_position > 0.6:
                            wheat_bonus += 10
                    
                    total_score = base_score + tm_bonus + dimer_penalty + size_bonus + tm_range_bonus + wheat_bonus
                    total_score = max(0, min(100, total_score))
                    
                    # 最终可用性检查
                    is_usable = (
                        fwd_eval['score'] >= 50 and
                        rev_eval['score'] >= 50 and
                        tm_diff <= config.MAX_TM_DIFF and
                        not has_dimer and
                        config.MIN_TM <= fwd_eval['tm'] <= config.MAX_TM and
                        config.MIN_TM <= rev_eval['tm'] <= config.MAX_TM
                    )
                    
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
                        'total_score': total_score,
                        'is_usable': is_usable,
                        'wheat_mode': config.WHEAT_MODE,
                        'wheat_issues': wheat_issues if config.WHEAT_MODE else [],
                        'fwd_position_percent': round(fwd_start / seq_len * 100, 1),
                        'rev_position_percent': round(rev_end / seq_len * 100, 1)
                    }
                    all_pairs.append(pair)
    
    # 排序：优先可用的，然后按评分
    all_pairs.sort(key=lambda x: (x['is_usable'], x['total_score']), reverse=True)
    
    # 去重并选择最佳方案
    unique_pairs = []
    seen = set()
    for pair in all_pairs:
        # 使用序列作为去重键
        key = (pair['forward'], pair['reverse'])
        if key not in seen:
            seen.add(key)
            unique_pairs.append(pair)
            if len(unique_pairs) >= num_pairs * 3:  # 先收集更多候选
                break
    
    # 确保选择的引物对有足够的位置多样性
    final_pairs = []
    used_positions = []
    
    for pair in unique_pairs:
        # 检查位置是否与已选引物太接近
        fwd_pos = pair['fwd_start']
        too_close = False
        for pos in used_positions:
            if abs(fwd_pos - pos) < 15:  # 至少间隔15bp
                too_close = True
                break
        
        if not too_close:
            final_pairs.append(pair)
            used_positions.append(fwd_pos)
            if len(final_pairs) >= num_pairs:
                break
    
    # 如果位置多样性导致数量不足，补充更多
    if len(final_pairs) < num_pairs:
        for pair in unique_pairs:
            if pair not in final_pairs:
                final_pairs.append(pair)
                if len(final_pairs) >= num_pairs:
                    break
    
    return final_pairs


# ==================== CSV导出函数 ====================

def generate_kasp_csv(schemes: List[Dict], seq_id: str) -> str:
    """生成KASP引物CSV内容"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    writer.writerow(['KASP引物设计报告'])
    writer.writerow(['序列ID', seq_id])
    writer.writerow(['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
    
    # 检查是否为小麦模式
    is_wheat_mode = schemes[0].get('wheat_mode', False) if schemes else False
    if is_wheat_mode:
        writer.writerow(['模式', '🌾 小麦KASP模式 (五大忌检测)'])
    writer.writerow([])
    
    # 表头
    headers = ['方案', '评分', '等级', '可用性',
               'FAM引物(完整)', 'HEX引物(完整)', '通用反向引物',
               'Allele1-Tm', 'Allele2-Tm', 'Tm差异',
               'Allele1-GC%', 'Allele2-GC%', 'Rev-Tm', 'Rev-GC%',
               '产物大小', '错配位置', '错配变化']
    
    if is_wheat_mode:
        headers.extend(['小麦评估', '注意事项'])
    
    writer.writerow(headers)
    
    for i, scheme in enumerate(schemes, 1):
        grade, stars, _ = get_quality_grade(scheme['total_score'])
        is_usable = scheme.get('is_usable', True)
        
        row = [
            f"方案{i}", f"{scheme['total_score']:.1f}", f"{grade} {stars}",
            "推荐" if is_usable else "慎用",
            scheme['fwd_allele1_full'], scheme['fwd_allele2_full'], scheme['reverse'],
            f"{scheme['eval_fwd1']['tm']}°C", f"{scheme['eval_fwd2']['tm']}°C", f"{scheme['tm_diff']:.1f}°C",
            f"{scheme['eval_fwd1']['gc_content']:.1f}%", f"{scheme['eval_fwd2']['gc_content']:.1f}%",
            f"{scheme['eval_rev']['tm']}°C", f"{scheme['eval_rev']['gc_content']:.1f}%",
            f"{scheme['product_size']}bp", f"n{scheme['mismatch_pos']}", scheme['mismatch_change']
        ]
        
        if is_wheat_mode:
            amplicon_status = scheme.get('wheat_details', {}).get('amplicon_status', '')
            wheat_issues = scheme.get('wheat_issues', [])
            row.append(amplicon_status)
            row.append("; ".join(wheat_issues) if wheat_issues else "无")
        
        writer.writerow(row)
    
    # 小麦模式提醒
    if is_wheat_mode:
        writer.writerow([])
        writer.writerow(['=== 小麦KASP五大忌提醒 ==='])
        writer.writerow(['1', '大忌#1 同源基因干扰', '请将引物BLAST到A/B/D三个基因组验证特异性'])
        writer.writerow(['2', '大忌#2 侧翼SNP干扰', '确认Allele引物区域无其他SNP/InDel'])
        writer.writerow(['3', '大忌#3 扩增子过长', '已优化为50-100bp'])
        writer.writerow(['4', '大忌#4 GC极端', '已检测30-65%范围'])
        writer.writerow(['5', '大忌#5 重复序列', '已检测转座子/SSR'])
        writer.writerow(['推荐工具', 'PolyMarker', 'http://polymarker.tgac.ac.uk/'])
        writer.writerow(['推荐工具', 'CerealsDB', 'http://www.cerealsdb.uk.net/'])
    
    return output.getvalue()


def generate_regular_csv(pairs: List[Dict], seq_id: str) -> str:
    """生成常规PCR引物CSV内容"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    writer.writerow(['常规PCR引物设计报告'])
    writer.writerow(['序列ID', seq_id])
    writer.writerow(['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
    writer.writerow([])
    
    writer.writerow(['引物对', '评分', '等级', '可用性',
                     '正向引物', '反向引物',
                     'Fwd-Tm', 'Rev-Tm', 'Tm差异',
                     'Fwd-GC%', 'Rev-GC%',
                     'Fwd评分', 'Rev评分',
                     'Fwd位置', 'Rev位置', '产物大小', '二聚体风险', '注意事项'])
    
    for i, pair in enumerate(pairs, 1):
        grade, stars, _ = get_quality_grade(pair['total_score'])
        is_usable = pair.get('is_usable', True)
        all_issues = pair['fwd_eval']['issues'] + pair['rev_eval']['issues']
        if pair['has_dimer']:
            all_issues.append("引物二聚体")
        
        writer.writerow([
            f"引物对{i}", f"{pair['total_score']:.1f}", f"{grade} {stars}",
            "推荐" if is_usable else "慎用",
            pair['forward'], pair['reverse'],
            f"{pair['fwd_eval']['tm']}°C", f"{pair['rev_eval']['tm']}°C", f"{pair['tm_diff']:.1f}°C",
            f"{pair['fwd_eval']['gc_content']:.1f}%", f"{pair['rev_eval']['gc_content']:.1f}%",
            f"{pair['fwd_eval']['score']:.0f}", f"{pair['rev_eval']['score']:.0f}",
            f"{pair['fwd_start']}-{pair['fwd_end']}", f"{pair['rev_start']}-{pair['rev_end']}",
            f"{pair['product_size']}bp",
            "有" if pair['has_dimer'] else "无",
            "; ".join(set(all_issues)) if all_issues else "无"
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
        
        # 🌾 小麦模式
        st.markdown("---")
        wheat_mode = st.checkbox("🌾 小麦KASP模式", value=False,
                                  help="针对小麦六倍体(AABBDD)优化，检测五大忌")
        
        if wheat_mode:
            st.warning("""**🌾 小麦五大忌检测已启用：**
1️⃣ 同源基因干扰 → 请BLAST验证
2️⃣ 侧翼SNP干扰 → 自动检测
3️⃣ 扩增子过长 → 限制50-100bp
4️⃣ GC含量极端 → 检测30-65%
5️⃣ 重复序列 → 检测转座子/SSR""")
        
        with st.expander("高级参数"):
            min_primer_len = st.number_input("最小引物长度", 15, 25, 18)
            max_primer_len = st.number_input("最大引物长度", 20, 35, 30)
            min_tm = st.number_input("最低Tm (°C)", 50.0, 60.0, 55.0)
            max_tm = st.number_input("最高Tm (°C)", 60.0, 75.0, 68.0)
            
            if wheat_mode:
                st.markdown("**小麦专用参数**")
                product_min = st.number_input("最小产物(bp)", 40, 80, 50, help="KASP建议50-100bp")
                product_max = st.number_input("最大产物(bp)", 80, 150, 100, help="小麦建议≤100bp")
            else:
                product_min = 80
                product_max = 200
    
    if st.button("🧬 设计KASP引物", type="primary", use_container_width=True):
        try:
            with st.spinner("正在设计引物..."):
                upstream, downstream, allele1, allele2 = parse_snp_sequence(seq_input)
                
                config = KASPConfig()
                config.MIN_PRIMER_LEN = min_primer_len
                config.MAX_PRIMER_LEN = max_primer_len
                config.MIN_TM = min_tm
                config.MAX_TM = max_tm
                
                # 小麦模式配置
                config.WHEAT_MODE = wheat_mode
                if wheat_mode:
                    config.PRODUCT_MIN = product_min
                    config.PRODUCT_MAX = product_max
                    config.REV_MIN_DISTANCE = 30
                    config.REV_MAX_DISTANCE = 80
                
                schemes = design_kasp_primers_multi(upstream, downstream, allele1, allele2, config, num_schemes)
            
            if not schemes:
                st.error("❌ 未能设计出合适的引物，请检查序列或调整参数")
                return
            
            st.success(f"✅ 成功设计 {len(schemes)} 套引物方案！")
            
            # 小麦模式警告
            if wheat_mode:
                st.info("""**🌾 重要提醒 - 大忌#1 同源基因干扰：**
请将设计的引物序列BLAST到小麦A、B、D三个基因组，确认：
- 引物是否只在目标基因组有完美匹配
- 是否需要在Common引物区域使用Genome-specific SNP
- 推荐工具：[PolyMarker](http://polymarker.tgac.ac.uk/) | [CerealsDB](http://www.cerealsdb.uk.net/)""")
            
            # 显示SNP信息
            st.markdown(f"""
            **SNP信息：** `{allele1}/{allele2}` | 
            **上游序列：** {len(upstream)}bp | 
            **下游序列：** {len(downstream)}bp
            """)
            
            # 显示每个方案
            for i, scheme in enumerate(schemes, 1):
                grade, stars, css_class = get_quality_grade(scheme['total_score'])
                is_usable = scheme.get('is_usable', True)
                usable_icon = "✅" if is_usable else "⚠️"
                
                with st.expander(f"{usable_icon} 方案 #{i} | 评分: {scheme['total_score']:.1f} | {grade} {stars}", expanded=(i==1)):
                    # 可用性提示
                    if wheat_mode:
                        if is_usable:
                            st.success("✅ 该方案通过小麦五大忌检测，推荐使用")
                        else:
                            st.error("⚠️ 该方案存在小麦特异性问题，请谨慎使用或选择其他方案")
                    
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
                        
                        # 产物大小评估（小麦模式）
                        if wheat_mode:
                            amplicon_status = scheme.get('wheat_details', {}).get('amplicon_status', '')
                            st.write(f"产物大小: {scheme['product_size']} bp {amplicon_status}")
                        else:
                            st.write(f"产物大小: {scheme['product_size']} bp")
                    
                    # 显示问题
                    all_issues = scheme['eval_fwd1']['issues'] + scheme['eval_fwd2']['issues'] + scheme['eval_rev']['issues']
                    if scheme['has_dimer']:
                        all_issues.append("引物间可能形成二聚体")
                    
                    # 小麦特异性问题
                    wheat_issues = scheme.get('wheat_issues', [])
                    if wheat_issues:
                        all_issues.extend(wheat_issues)
                    
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
        
        # 小麦模式
        st.markdown("---")
        wheat_mode = st.checkbox("🌾 小麦特异性模式", value=False, 
                                  help="针对小麦A/B/D同源基因优化，避开5'端保守区")
        
        if wheat_mode:
            st.info("""**小麦模式已启用：**
- 自动避开5'端保守区域
- 优先在3'端/UTR区设计引物
- 评估同源基因特异性""")
            
            with st.expander("🌾 小麦参数设置"):
                avoid_5prime = st.slider("避开5'端区域(%)", 20, 60, 40,
                                         help="避开序列5'端的百分比，该区域同源基因通常高度保守")
                prefer_3prime = st.checkbox("优先3'端区域", value=True,
                                            help="3'UTR区域通常变异更多，有利于特异性扩增")
        else:
            avoid_5prime = 40
            prefer_3prime = True
        
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
                
                # 小麦模式参数
                config.WHEAT_MODE = wheat_mode
                if wheat_mode:
                    config.WHEAT_AVOID_5PRIME_PERCENT = avoid_5prime
                    config.WHEAT_PREFER_3PRIME = prefer_3prime
                
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
                is_usable = pair.get('is_usable', True)
                usable_icon = "✅" if is_usable else "⚠️"
                
                with st.expander(f"{usable_icon} 引物对 #{i} | 评分: {pair['total_score']:.1f} | {grade} {stars}", expanded=(i==1)):
                    # 可用性提示
                    if is_usable:
                        st.success("✅ 该引物对通过所有质量检测，推荐使用")
                    else:
                        st.warning("⚠️ 该引物对存在一些问题，建议优先选择其他方案")
                    
                    col_a, col_b = st.columns(2)
                    
                    with col_a:
                        st.markdown("**正向引物 (Forward)**")
                        st.code(f"5'- {pair['forward']} -3'")
                        fwd_score_grade, _, _ = get_quality_grade(pair['fwd_eval']['score'])
                        st.caption(f"位置: {pair['fwd_start']}-{pair['fwd_end']} | {len(pair['forward'])}bp | Tm: {pair['fwd_eval']['tm']}°C | GC: {pair['fwd_eval']['gc_content']:.1f}% | 评分: {pair['fwd_eval']['score']:.0f}({fwd_score_grade})")
                    
                    with col_b:
                        st.markdown("**反向引物 (Reverse)**")
                        st.code(f"5'- {pair['reverse']} -3'")
                        rev_score_grade, _, _ = get_quality_grade(pair['rev_eval']['score'])
                        st.caption(f"位置: {pair['rev_start']}-{pair['rev_end']} | {len(pair['reverse'])}bp | Tm: {pair['rev_eval']['tm']}°C | GC: {pair['rev_eval']['gc_content']:.1f}% | 评分: {pair['rev_eval']['score']:.0f}({rev_score_grade})")
                    
                    tm_status = "✓ 优秀" if pair['tm_diff'] <= 1.0 else ("△ 可接受" if pair['tm_diff'] <= 2.0 else "✗ 偏大")
                    dimer_status = "✗ 有风险" if pair['has_dimer'] else "✓ 无风险"
                    st.write(f"**Tm差异:** {pair['tm_diff']:.1f}°C ({tm_status}) | **产物大小:** {pair['product_size']} bp | **二聚体:** {dimer_status}")
                    
                    # 小麦模式：显示位置信息
                    if pair.get('wheat_mode'):
                        fwd_pos = pair.get('fwd_position_percent', 0)
                        rev_pos = pair.get('rev_position_percent', 0)
                        
                        # 位置评估
                        if fwd_pos < 30:
                            pos_status = "⚠️ 5'端保守区"
                        elif fwd_pos < 50:
                            pos_status = "△ 中间区域"
                        else:
                            pos_status = "✓ 3'端变异区"
                        
                        st.write(f"**🌾 位置:** Fwd在{fwd_pos:.0f}% | Rev在{rev_pos:.0f}% ({pos_status})")
                    
                    all_issues = pair['fwd_eval']['issues'] + pair['rev_eval']['issues']
                    if pair['has_dimer']:
                        all_issues.append("引物间可能形成二聚体")
                    
                    # 小麦特异性问题
                    wheat_issues = pair.get('wheat_issues', [])
                    if wheat_issues:
                        all_issues.extend(wheat_issues)
                    
                    if all_issues:
                        st.warning("⚠️ 注意事项: " + " | ".join(set(all_issues)))
            
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
    """引物分析工具 - 增强版"""
    st.markdown("### 🔍 引物质量分析")
    
    st.markdown("""
    <div class="info-box">
    <b>功能说明：</b><br>
    - 支持单条引物分析<br>
    - 支持引物对分析（两条引物间的相互作用）<br>
    - 支持小麦KASP引物特异性分析<br>
    - 提供详细的质量评估和优化建议
    </div>
    """, unsafe_allow_html=True)
    
    # 选择分析模式
    analysis_mode = st.radio(
        "选择分析模式",
        ["单引物分析", "引物对分析", "小麦KASP引物分析"],
        horizontal=True
    )
    
    if analysis_mode == "单引物分析":
        st.markdown("#### 📝 输入引物序列")
        primer_input = st.text_area(
            "引物序列 (5'→3')",
            placeholder="例如: ATGCGATCGATCGATCGATCG\n只输入ATGC碱基，自动过滤其他字符",
            height=100
        )
        
        primer_name = st.text_input("引物名称（可选）", value="My_Primer")
        
        # 分析类型选择
        col1, col2 = st.columns(2)
        with col1:
            primer_type = st.selectbox("引物类型", ["常规PCR", "KASP Allele", "KASP Common"])
        with col2:
            check_wheat = st.checkbox("小麦特异性检测", value=False)
        
        if st.button("🔍 开始分析", type="primary"):
            if not primer_input:
                st.warning("请输入引物序列")
                return
            
            primer = re.sub(r'[^ATGC]', '', primer_input.upper())
            
            if len(primer) < 10:
                st.error("❌ 引物序列过短（<10bp），无法进行有效分析")
                return
            
            # 根据类型选择配置
            if primer_type == "KASP Allele" or primer_type == "KASP Common":
                config = KASPConfig()
                config.WHEAT_MODE = check_wheat
            else:
                config = RegularPCRConfig()
                config.WHEAT_MODE = check_wheat
            
            # 基础分析
            result = evaluate_primer_quality(primer, config)
            grade, stars, css_class = get_quality_grade(result['score'])
            
            # 显示结果
            st.markdown("---")
            st.markdown(f"### 📊 分析报告: {primer_name}")
            st.markdown(f"<h3 style='text-align:center;'><span class='{css_class}'>{grade} {stars} ({result['score']:.1f}分)</span></h3>", unsafe_allow_html=True)
            
            # 引物信息
            st.code(f"5'- {primer} -3'")
            st.caption(f"长度: {len(primer)} bp")
            
            # 指标展示
            col1, col2, col3, col4, col5 = st.columns(5)
            
            with col1:
                tm_color = "🟢" if 58 <= result['tm'] <= 62 else ("🟡" if 55 <= result['tm'] <= 68 else "🔴")
                st.metric("Tm值", f"{result['tm']}°C", delta=f"{tm_color}")
            
            with col2:
                gc = result['gc_content']
                gc_color = "🟢" if 40 <= gc <= 60 else ("🟡" if 30 <= gc <= 70 else "🔴")
                st.metric("GC含量", f"{gc:.1f}%", delta=f"{gc_color}")
            
            with col3:
                len_color = "🟢" if 18 <= len(primer) <= 25 else ("🟡" if 15 <= len(primer) <= 30 else "🔴")
                st.metric("长度", f"{len(primer)}bp", delta=f"{len_color}")
            
            with col4:
                end_base = primer[-1]
                end_color = "🟢" if end_base in ['G', 'C'] else "🟡"
                st.metric("3'端", end_base, delta=f"{end_color}")
            
            with col5:
                st.metric("综合评分", f"{result['score']:.0f}", delta="/ 100")
            
            # 详细检测结果
            st.markdown("---")
            col_a, col_b = st.columns(2)
            
            with col_a:
                st.markdown("**🔬 结构检测**")
                hairpin_status = "❌ 检测到" if result['has_hairpin'] else "✅ 未检测到"
                st.write(f"**发夹结构:** {hairpin_status}")
                if result['has_hairpin']:
                    st.caption("⚠️ 可能影响引物特异性和扩增效率")
                
                dimer_status = "❌ 有风险" if result['has_self_dimer'] else "✅ 无风险"
                st.write(f"**自身二聚体:** {dimer_status}")
                if result['has_self_dimer']:
                    st.caption("⚠️ 可能导致引物-引物扩增")
                
                # 检查重复
                has_repeat = check_repeat_region(primer)
                repeat_status = "❌ 检测到" if has_repeat else "✅ 无"
                st.write(f"**重复序列:** {repeat_status}")
            
            with col_b:
                st.markdown("**🎯 3'端分析**")
                three_prime_icon = "✅" if result['three_prime_ok'] else "⚠️"
                st.write(f"**3'端稳定性:** {three_prime_icon} {result['three_prime_msg']}")
                
                end_5 = primer[-5:]
                gc_end = end_5.count('G') + end_5.count('C')
                st.write(f"**3'端5bp:** `{end_5}` (GC={gc_end}/5)")
                
                if gc_end > 3:
                    st.caption("⚠️ 3'端GC过多，可能非特异性结合")
                elif gc_end < 1:
                    st.caption("⚠️ 3'端GC过少，结合不稳定")
                else:
                    st.caption("✓ 3'端GC含量良好")
            
            # 小麦特异性分析
            if check_wheat:
                st.markdown("---")
                st.markdown("**🌾 小麦特异性分析**")
                
                # 检测重复序列
                has_wheat_repeat, wheat_repeat_issues = check_wheat_repeat_sequences(primer)
                
                if has_wheat_repeat:
                    st.error("❌ 检测到小麦常见重复序列特征")
                    for issue in wheat_repeat_issues:
                        st.write(f"  - {issue}")
                else:
                    st.success("✅ 未检测到明显的重复序列特征")
                
                # GC检测
                gc_extreme, gc_val, gc_msg = check_gc_extreme(primer, strict=True)
                if gc_extreme:
                    st.error(f"❌ {gc_msg}")
                else:
                    st.success(f"✅ GC含量符合小麦KASP标准")
                
                # 序列复杂度
                complexity = analyze_sequence_complexity(primer)
                complexity_score = complexity['complexity_score']
                
                if complexity_score >= 75:
                    st.success(f"✅ 序列复杂度高 ({complexity_score:.0f}/100) - 特异性好")
                elif complexity_score >= 60:
                    st.info(f"ℹ️ 序列复杂度中等 ({complexity_score:.0f}/100)")
                else:
                    st.warning(f"⚠️ 序列复杂度低 ({complexity_score:.0f}/100) - 可能多位点匹配")
            
            # 问题汇总
            if result['issues']:
                st.markdown("---")
                st.markdown("**⚠️ 检测到的问题**")
                for issue in result['issues']:
                    st.write(f"• {issue}")
            
            # 优化建议
            st.markdown("---")
            st.markdown("**💡 优化建议**")
            suggestions = []
            
            if result['tm'] < 55:
                suggestions.append("增加引物长度或提高GC含量以提升Tm值")
            elif result['tm'] > 68:
                suggestions.append("缩短引物或降低GC含量以降低Tm值")
            
            if result['gc_content'] < 40:
                suggestions.append("选择GC含量更高的区域设计引物")
            elif result['gc_content'] > 60:
                suggestions.append("选择GC含量更适中的区域设计引物")
            
            if result['has_hairpin']:
                suggestions.append("改变引物位置或长度以避免发夹结构")
            
            if result['has_self_dimer']:
                suggestions.append("调整引物序列以减少自身互补区域")
            
            if not result['three_prime_ok']:
                suggestions.append("调整引物使3'端有1-2个G或C")
            
            if check_wheat and has_wheat_repeat:
                suggestions.append("更换设计区域，避开重复序列/转座子区域")
            
            if suggestions:
                for i, sug in enumerate(suggestions, 1):
                    st.write(f"{i}. {sug}")
            else:
                st.success("✅ 引物质量良好，无需优化")
    
    elif analysis_mode == "引物对分析":
        st.markdown("#### 📝 输入引物对序列")
        
        col1, col2 = st.columns(2)
        with col1:
            fwd_input = st.text_area("正向引物 (Forward 5'→3')", 
                                     placeholder="ATGCGATCGATCGATCG",
                                     height=100)
            fwd_name = st.text_input("正向引物名称", value="Forward_Primer")
        
        with col2:
            rev_input = st.text_area("反向引物 (Reverse 5'→3')", 
                                     placeholder="CGATCGATCGATCGAT",
                                     height=100)
            rev_name = st.text_input("反向引物名称", value="Reverse_Primer")
        
        if st.button("🔍 分析引物对", type="primary"):
            if not fwd_input or not rev_input:
                st.warning("请输入两条引物序列")
                return
            
            fwd = re.sub(r'[^ATGC]', '', fwd_input.upper())
            rev = re.sub(r'[^ATGC]', '', rev_input.upper())
            
            if len(fwd) < 10 or len(rev) < 10:
                st.error("❌ 引物序列过短")
                return
            
            # 分析两条引物
            config = RegularPCRConfig()
            fwd_result = evaluate_primer_quality_strict(fwd, config)
            rev_result = evaluate_primer_quality_strict(rev, config)
            
            # 引物对分析
            tm_diff = abs(fwd_result['tm'] - rev_result['tm'])
            has_dimer = check_primer_dimer(fwd, rev)
            
            # 综合评分
            pair_score = (fwd_result['score'] + rev_result['score']) / 2
            if tm_diff <= 1.0:
                pair_score += 10
            elif tm_diff <= 2.0:
                pair_score += 5
            else:
                pair_score -= 10
            
            if has_dimer:
                pair_score -= 15
            
            pair_score = max(0, min(100, pair_score))
            grade, stars, css_class = get_quality_grade(pair_score)
            
            # 显示结果
            st.markdown("---")
            st.markdown(f"### 📊 引物对分析报告")
            st.markdown(f"<h3 style='text-align:center;'><span class='{css_class}'>综合评分: {grade} {stars} ({pair_score:.1f}分)</span></h3>", unsafe_allow_html=True)
            
            # 两条引物对比
            col_a, col_b = st.columns(2)
            
            with col_a:
                st.markdown(f"**正向引物: {fwd_name}**")
                st.code(f"5'- {fwd} -3'")
                fwd_grade, fwd_stars, _ = get_quality_grade(fwd_result['score'])
                st.caption(f"评分: {fwd_result['score']:.0f} ({fwd_grade} {fwd_stars})")
                
                st.write(f"**长度:** {len(fwd)} bp")
                st.write(f"**Tm:** {fwd_result['tm']}°C")
                st.write(f"**GC:** {fwd_result['gc_content']:.1f}%")
                
                if fwd_result['issues']:
                    st.warning("问题: " + ", ".join(fwd_result['issues'][:2]))
            
            with col_b:
                st.markdown(f"**反向引物: {rev_name}**")
                st.code(f"5'- {rev} -3'")
                rev_grade, rev_stars, _ = get_quality_grade(rev_result['score'])
                st.caption(f"评分: {rev_result['score']:.0f} ({rev_grade} {rev_stars})")
                
                st.write(f"**长度:** {len(rev)} bp")
                st.write(f"**Tm:** {rev_result['tm']}°C")
                st.write(f"**GC:** {rev_result['gc_content']:.1f}%")
                
                if rev_result['issues']:
                    st.warning("问题: " + ", ".join(rev_result['issues'][:2]))
            
            # 配对分析
            st.markdown("---")
            st.markdown("**🔗 配对分析**")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                tm_status = "✅ 优秀" if tm_diff <= 1.0 else ("⚠️ 可接受" if tm_diff <= 2.0 else "❌ 过大")
                st.metric("Tm差异", f"{tm_diff:.1f}°C", delta=tm_status)
            
            with col2:
                dimer_status = "❌ 检测到" if has_dimer else "✅ 无风险"
                st.metric("引物二聚体", dimer_status)
            
            with col3:
                avg_tm = (fwd_result['tm'] + rev_result['tm']) / 2
                st.metric("平均Tm", f"{avg_tm:.1f}°C")
            
            # 详细说明
            if has_dimer:
                st.error("""**❌ 检测到引物二聚体风险**
- 两条引物之间可能形成互补配对
- 可能导致引物-引物扩增而非目标扩增
- 建议：更改其中一条引物的序列""")
            else:
                st.success("✅ 未检测到明显的引物二聚体风险")
            
            if tm_diff > 2.0:
                st.warning(f"""**⚠️ Tm差异过大 ({tm_diff:.1f}°C)**
- 两条引物退火温度相差较大
- 可能导致扩增效率不平衡
- 建议：调整引物使Tm差异 ≤2°C""")
            
            # 建议退火温度
            st.markdown("---")
            st.markdown("**🌡️ 建议PCR条件**")
            recommended_tm = min(fwd_result['tm'], rev_result['tm']) - 5
            st.write(f"**推荐退火温度:** {recommended_tm:.0f}°C (较低Tm - 5°C)")
            st.write(f"**梯度PCR范围:** {recommended_tm-3:.0f}°C ~ {recommended_tm+3:.0f}°C")
    
    else:  # 小麦KASP引物分析
        st.markdown("#### 📝 输入KASP引物组")
        
        st.info("""**小麦KASP引物组包括：**
- 2条等位基因特异性引物（带FAM/HEX荧光尾巴）
- 1条通用反向引物（Common Primer）""")
        
        allele1_input = st.text_area("Allele 1 引物（完整，含FAM尾巴）", height=80)
        allele2_input = st.text_area("Allele 2 引物（完整，含HEX尾巴）", height=80)
        common_input = st.text_area("Common 反向引物", height=80)
        
        if st.button("🌾 分析小麦KASP引物", type="primary"):
            if not (allele1_input and allele2_input and common_input):
                st.warning("请输入完整的KASP引物组")
                return
            
            allele1 = re.sub(r'[^ATGC]', '', allele1_input.upper())
            allele2 = re.sub(r'[^ATGC]', '', allele2_input.upper())
            common = re.sub(r'[^ATGC]', '', common_input.upper())
            
            config = KASPConfig()
            config.WHEAT_MODE = True
            
            # 提取核心序列（去除荧光尾巴）
            fam_tail_len = len(config.FAM_TAIL)
            hex_tail_len = len(config.HEX_TAIL)
            
            core1 = allele1[fam_tail_len:] if len(allele1) > fam_tail_len else allele1
            core2 = allele2[hex_tail_len:] if len(allele2) > hex_tail_len else allele2
            
            # 分析核心序列
            eval1 = evaluate_primer_quality(core1, config)
            eval2 = evaluate_primer_quality(core2, config)
            eval_common = evaluate_primer_quality(common, config)
            
            # 小麦特异性分析
            upstream = ""  # 简化，仅分析引物本身
            downstream = ""
            wheat_bonus, wheat_issues, wheat_details = evaluate_kasp_wheat_specificity(
                upstream, downstream, core1, common, config
            )
            
            # 显示结果
            st.markdown("---")
            st.markdown("### 🌾 小麦KASP引物分析报告")
            
            # 五大忌检测结果
            st.markdown("**五大忌检测结果：**")
            checks = [
                ("1️⃣ 同源基因干扰", "需要BLAST验证", "⚠️"),
                ("2️⃣ 侧翼SNP干扰", wheat_details.get('flanking_risk', False), "❌" if wheat_details.get('flanking_risk', False) else "✅"),
                ("3️⃣ 扩增子长度", "无法计算（缺少完整序列）", "ℹ️"),
                ("4️⃣ GC含量极端", f"Core1:{eval1['gc_content']:.1f}% Core2:{eval2['gc_content']:.1f}% Common:{eval_common['gc_content']:.1f}%", 
                 "✅" if all(30 <= x['gc_content'] <= 65 for x in [eval1, eval2, eval_common]) else "❌"),
                ("5️⃣ 重复序列", wheat_details.get('has_repeat', False), "❌" if wheat_details.get('has_repeat', False) else "✅"),
            ]
            
            for check_name, check_result, check_icon in checks:
                if isinstance(check_result, bool):
                    result_text = "检测到" if check_result else "通过"
                else:
                    result_text = str(check_result)
                st.write(f"{check_icon} **{check_name}:** {result_text}")
            
            # 引物详情
            st.markdown("---")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("**Allele 1 引物**")
                st.caption(f"核心: {len(core1)}bp | Tm: {eval1['tm']}°C | GC: {eval1['gc_content']:.1f}%")
                grade1, _, _ = get_quality_grade(eval1['score'])
                st.write(f"评分: {eval1['score']:.0f} ({grade1})")
            
            with col2:
                st.markdown("**Allele 2 引物**")
                st.caption(f"核心: {len(core2)}bp | Tm: {eval2['tm']}°C | GC: {eval2['gc_content']:.1f}%")
                grade2, _, _ = get_quality_grade(eval2['score'])
                st.write(f"评分: {eval2['score']:.0f} ({grade2})")
            
            with col3:
                st.markdown("**Common 引物**")
                st.caption(f"长度: {len(common)}bp | Tm: {eval_common['tm']}°C | GC: {eval_common['gc_content']:.1f}%")
                grade_c, _, _ = get_quality_grade(eval_common['score'])
                st.write(f"评分: {eval_common['score']:.0f} ({grade_c})")
            
            # Tm匹配
            tm_diff = abs(eval1['tm'] - eval2['tm'])
            st.markdown("---")
            if tm_diff <= 1.0:
                st.success(f"✅ Allele引物Tm差异: {tm_diff:.1f}°C (优秀)")
            elif tm_diff <= 2.0:
                st.info(f"ℹ️ Allele引物Tm差异: {tm_diff:.1f}°C (可接受)")
            else:
                st.warning(f"⚠️ Allele引物Tm差异: {tm_diff:.1f}°C (过大)")
            
            # 小麦特异性问题
            if wheat_issues:
                st.markdown("---")
                st.markdown("**⚠️ 小麦特异性问题：**")
                for issue in wheat_issues:
                    st.write(f"• {issue}")
            
            # 建议
            st.markdown("---")
            st.markdown("**💡 重要建议：**")
            st.write("1. 将所有引物序列BLAST到小麦A、B、D三个基因组")
            st.write("2. 确认引物只匹配目标基因组")
            st.write("3. 如需基因组特异性，在Common引物区域添加Homoeologous SNP")
            st.write("4. 使用PolyMarker工具验证设计")
            st.write("5. 推荐产物大小50-100bp")


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
