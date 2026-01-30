"""
KASP & 常规PCR 引物设计工具 - Streamlit Web版
版本: v7.0 Web (Primer3-py重构版)
功能: KASP引物设计、常规PCR引物设计、质量评估、CSV导出
重构: 使用primer3-py库进行专业的热力学计算，参考polyoligo-kasp设计理念
"""

import streamlit as st
import re
import csv
import io
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass

# 尝试导入primer3库
try:
    import primer3
    PRIMER3_AVAILABLE = True
    PRIMER3_VERSION = getattr(primer3, '__version__', 'unknown')
except ImportError:
    PRIMER3_AVAILABLE = False
    PRIMER3_VERSION = None

# ==================== 页面配置 ====================
st.set_page_config(
    page_title="引物设计工具 v7.0",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 显示Primer3状态 (在侧边栏底部)
def show_primer3_status():
    """显示Primer3库状态"""
    if PRIMER3_AVAILABLE:
        st.sidebar.success(f"✅ Primer3-py v{PRIMER3_VERSION} 已加载")
    else:
        st.sidebar.error("""⚠️ primer3-py未安装
        
请运行以下命令安装：
```bash
pip install primer3-py
```
当前使用内置算法（精度较低）""")


def show_kasp_features_info():
    """显示KASP引物设计的核心功能说明"""
    with st.sidebar.expander("🔬 KASP设计核心功能", expanded=False):
        st.markdown("""
        ### 1️⃣ 人工错配 (Deliberate Mismatch)
        
        **LGC标准** - 在ASP引物3'端倒数第3位(n-3)引入人工错配
        
        - **强SNP (G/C)**: 使用强destabilizing错配
        - **弱SNP (A/T)**: 使用中等destabilizing错配
        
        **作用**: 增强等位基因特异性，减少非特异性扩增
        
        ---
        
        ### 2️⃣ 救援模式 (Rescue Mode)
        
        自动处理**AT-rich序列**的两轮筛选机制:
        
        | 参数 | 标准模式 | 救援模式 |
        |------|---------|---------|
        | 最大长度 | 25bp | 32bp |
        | GC下限 | 35% | 20% |
        | Tm下限 | 55°C | 52°C |
        
        **触发**: 标准参数无法设计时自动启用
        
        ---
        
        ### 3️⃣ Tm平衡优化
        
        - **ASP引物**: 只计算核心序列Tm(不含FAM/HEX)
        - **Common引物**: 目标Tm 60-62°C
        - **匹配策略**: Common引物长度可达30bp以匹配ASP的Tm
        
        ---
        
        ### 4️⃣ 智能失败诊断
        
        当救援模式也失败时，自动分析原因：
        
        | 诊断项 | 阈值 |
        |--------|------|
        | GC极低 | <25% → 建议换SNP |
        | GC偏低 | <35% → 建议加长引物 |
        | 发夹风险 | Tm>45°C → 避开区域 |
        | 重复序列 | SSR/转座子 → 警告 |
        | 序列过短 | <25bp → 提供更长序列 |
        """)


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
    
    # === 标准模式参数 ===
    MIN_PRIMER_LEN: int = 18
    MAX_PRIMER_LEN: int = 25      # 标准模式最大长度
    OPTIMAL_PRIMER_LEN: int = 20
    MIN_TM: float = 55.0
    MAX_TM: float = 68.0
    OPTIMAL_TM: float = 62.0
    MIN_GC: float = 35.0          # 标准模式GC下限
    MAX_GC: float = 65.0          # 小麦大忌：避免>65%
    OPTIMAL_GC_MIN: float = 40.0
    OPTIMAL_GC_MAX: float = 55.0
    MAX_TM_DIFF: float = 2.0
    
    # KASP产物大小：小麦建议50-100bp
    REV_MIN_DISTANCE: int = 30   # 最近距离，产物约50bp
    REV_MAX_DISTANCE: int = 80   # 最远距离，产物约100bp
    PRODUCT_MIN: int = 50        # 小麦大忌#3：扩增子过长
    PRODUCT_MAX: int = 120       # KASP建议50-100bp
    MISMATCH_POSITIONS: List[int] = None
    
    # === 救援模式参数 (AT-rich序列) ===
    RESCUE_MODE_ENABLED: bool = True       # 是否启用救援模式
    RESCUE_MAX_PRIMER_LEN: int = 32        # 救援模式最大引物长度
    RESCUE_MIN_GC: float = 20.0            # 救援模式GC下限
    RESCUE_MIN_TM: float = 52.0            # 救援模式Tm下限
    
    # === Common引物Tm平衡参数 ===
    COMMON_TARGET_TM: float = 61.0         # Common引物目标Tm
    COMMON_TM_TOLERANCE: float = 1.5       # Common引物Tm容差
    COMMON_MAX_LEN: int = 30               # Common引物最大长度(用于Tm平衡)
    ASP_COMMON_TM_DIFF_MAX: float = 3.0    # ASP与Common的最大Tm差
    
    # === Tm贪婪延伸优化参数 (LGC商业设计标准) ===
    TM_GREEDY_EXTENSION: bool = True       # 是否启用Tm贪婪延伸
    ASP_TARGET_TM_MIN: float = 60.0        # ASP引物目标Tm下限
    ASP_TARGET_TM_MAX: float = 62.0        # ASP引物目标Tm上限 (LGC标准)
    TM_GREEDY_TOLERANCE: float = 1.0       # Tm达标容差
    
    # 🌾 小麦特异性参数
    WHEAT_MODE: bool = False
    WHEAT_CHECK_FLANKING_SNP: bool = True   # 大忌#2：检查侧翼干扰SNP
    WHEAT_CHECK_REPEAT: bool = True         # 大忌#5：检查重复序列
    WHEAT_STRICT_GC: bool = True            # 大忌#4：GC含量严格模式
    
    def __post_init__(self):
        if self.MISMATCH_POSITIONS is None:
            self.MISMATCH_POSITIONS = [-3, -2, -4]  # 优先n-3位置
    
    def get_rescue_config(self) -> 'KASPConfig':
        """
        获取救援模式配置 (用于AT-rich序列)
        放宽参数以确保能设计出引物
        """
        rescue = KASPConfig(
            FAM_TAIL=self.FAM_TAIL,
            HEX_TAIL=self.HEX_TAIL,
            MIN_PRIMER_LEN=self.MIN_PRIMER_LEN,
            MAX_PRIMER_LEN=self.RESCUE_MAX_PRIMER_LEN,  # 增加到32bp
            OPTIMAL_PRIMER_LEN=24,                       # 最优长度增加
            MIN_TM=self.RESCUE_MIN_TM,                   # 降低Tm下限
            MAX_TM=self.MAX_TM,
            OPTIMAL_TM=58.0,                             # 降低最优Tm
            MIN_GC=self.RESCUE_MIN_GC,                   # 降低GC下限到20%
            MAX_GC=self.MAX_GC,
            OPTIMAL_GC_MIN=30.0,
            OPTIMAL_GC_MAX=55.0,
            MAX_TM_DIFF=3.0,                             # 放宽Tm差异
            REV_MIN_DISTANCE=self.REV_MIN_DISTANCE,
            REV_MAX_DISTANCE=self.REV_MAX_DISTANCE + 20,  # 扩大搜索范围
            PRODUCT_MIN=self.PRODUCT_MIN,
            PRODUCT_MAX=150,                              # 允许更长产物
            MISMATCH_POSITIONS=self.MISMATCH_POSITIONS,
            RESCUE_MODE_ENABLED=False,                    # 防止递归
            WHEAT_MODE=self.WHEAT_MODE,
            WHEAT_CHECK_FLANKING_SNP=self.WHEAT_CHECK_FLANKING_SNP,
            WHEAT_CHECK_REPEAT=False,                     # 放宽重复检测
            WHEAT_STRICT_GC=False,                        # 放宽GC检测
        )
        return rescue


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


# ==================== Primer3 热力学参数配置 ====================

@dataclass
class Primer3ThermoParams:
    """Primer3热力学计算参数 - 参考polyoligo-kasp"""
    mv_conc: float = 50.0      # 单价阳离子浓度 (mM)
    dv_conc: float = 1.5       # 二价阳离子浓度 (mM)
    dntp_conc: float = 0.6     # dNTP浓度 (mM)
    dna_conc: float = 250.0    # DNA/引物浓度 (nM)
    temp_c: float = 37.0       # 模拟温度 (°C)
    max_loop: int = 30         # 最大环大小
    
    # KASP特定参数
    kasp_mv_conc: float = 50.0
    kasp_dv_conc: float = 2.0   # KASP通常使用更高的Mg2+
    kasp_dntp_conc: float = 0.8
    kasp_dna_conc: float = 200.0


# 全局热力学参数实例
THERMO_PARAMS = Primer3ThermoParams()


# ==================== 核心计算函数 (Primer3-py重构) ====================

def calc_gc_content(seq: str) -> float:
    """计算GC含量百分比"""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def calc_tm(seq: str, mv_conc: float = 50.0, dv_conc: float = 1.5, 
            dntp_conc: float = 0.6, dna_conc: float = 250.0) -> float:
    """
    使用Primer3计算Tm值 (SantaLucia方法)
    
    参数:
        seq: 引物序列
        mv_conc: 单价阳离子浓度 (mM)
        dv_conc: 二价阳离子浓度 (mM)  
        dntp_conc: dNTP浓度 (mM)
        dna_conc: DNA浓度 (nM)
    
    返回:
        Tm值 (°C)
    """
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)  # 清除非标准碱基
    
    if len(seq) < 5:
        # 序列过短，使用简单Wallace公式
        gc_count = seq.count('G') + seq.count('C')
        at_count = seq.count('A') + seq.count('T')
        return round(2 * at_count + 4 * gc_count, 1)
    
    if PRIMER3_AVAILABLE:
        try:
            tm = primer3.calc_tm(
                seq,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
                tm_method='santalucia',
                salt_corrections_method='santalucia'
            )
            return round(tm, 1)
        except Exception:
            pass
    
    # 回退到内置最近邻算法
    return _calc_tm_fallback(seq, mv_conc, dv_conc, dna_conc)


def _calc_tm_fallback(seq: str, na_conc: float = 50.0, 
                      mg_conc: float = 1.5, primer_conc: float = 250.0) -> float:
    """内置最近邻法Tm计算(回退方案)"""
    import math
    seq = seq.upper()
    length = len(seq)
    
    # SantaLucia 1998 统一参数
    nn_params = {
        'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
        'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
        'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
        'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
        'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
        'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
        'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
        'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
    }
    init_params = {
        'G': (0.1, -2.8), 'C': (0.1, -2.8),
        'A': (2.3, 4.1), 'T': (2.3, 4.1)
    }
    
    dH, dS = 0.0, 0.0
    if seq[0] in init_params:
        dH += init_params[seq[0]][0]
        dS += init_params[seq[0]][1]
    if seq[-1] in init_params:
        dH += init_params[seq[-1]][0]
        dS += init_params[seq[-1]][1]
    
    for i in range(length - 1):
        dinuc = seq[i:i+2]
        if dinuc in nn_params:
            dH += nn_params[dinuc][0]
            dS += nn_params[dinuc][1]
    
    na_equiv = na_conc + 120 * math.sqrt(mg_conc)
    na_molar = na_equiv / 1000.0
    dS_corrected = dS + 0.368 * (length - 1) * math.log(na_molar)
    
    R = 1.987
    ct = primer_conc * 1e-9
    tm = (dH * 1000) / (dS_corrected + R * math.log(ct / 4)) - 273.15
    
    return round(tm, 1)


# 保留旧函数名以兼容
calc_tm_nearest_neighbor = calc_tm
calc_tm_simple = lambda seq: calc_tm(seq)


def reverse_complement(seq: str) -> str:
    """生成反向互补序列"""
    if PRIMER3_AVAILABLE:
        try:
            from primer3 import p3helpers
            return p3helpers.reverse_complement(seq.upper())
        except Exception:
            pass
    
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


def check_hairpin(seq: str, min_stem: int = 4, min_loop: int = 3) -> Tuple[bool, Optional[float]]:
    """
    使用Primer3检测发夹结构
    
    返回: (是否有发夹, 发夹Tm值)
    """
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)
    
    if len(seq) < 10:
        return False, None
    
    if PRIMER3_AVAILABLE:
        try:
            result = primer3.calc_hairpin(
                seq,
                mv_conc=THERMO_PARAMS.mv_conc,
                dv_conc=THERMO_PARAMS.dv_conc,
                dntp_conc=THERMO_PARAMS.dntp_conc,
                dna_conc=THERMO_PARAMS.dna_conc,
                temp_c=THERMO_PARAMS.temp_c,
                max_loop=THERMO_PARAMS.max_loop
            )
            # 发夹Tm > 45°C 或 dG < -2 kcal/mol 认为有问题
            has_hairpin = result.structure_found and (result.tm > 45 or result.dg < -2000)
            return has_hairpin, result.tm if result.structure_found else None
        except Exception:
            pass
    
    # 回退到简单检测
    return _check_hairpin_fallback(seq, min_stem, min_loop), None


def _check_hairpin_fallback(seq: str, min_stem: int = 4, min_loop: int = 3) -> bool:
    """内置发夹检测(回退方案)"""
    seq = seq.upper()
    length = len(seq)
    
    for i in range(length - min_stem - min_loop):
        for stem_len in range(min_stem, min(8, (length - i - min_loop) // 2 + 1)):
            stem1 = seq[i:i + stem_len]
            for loop_len in range(min_loop, min(10, length - i - 2 * stem_len + 1)):
                j = i + stem_len + loop_len
                if j + stem_len <= length:
                    stem2 = seq[j:j + stem_len]
                    rc_stem2 = reverse_complement(stem2)[::-1]
                    if stem1 == rc_stem2:
                        return True
    return False


def check_homodimer(seq: str) -> Tuple[bool, Optional[float], Optional[float]]:
    """
    使用Primer3检测自身二聚体(同源二聚体)
    
    返回: (是否有二聚体风险, dG值, Tm值)
    """
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)
    
    if len(seq) < 10:
        return False, None, None
    
    if PRIMER3_AVAILABLE:
        try:
            result = primer3.calc_homodimer(
                seq,
                mv_conc=THERMO_PARAMS.mv_conc,
                dv_conc=THERMO_PARAMS.dv_conc,
                dntp_conc=THERMO_PARAMS.dntp_conc,
                dna_conc=THERMO_PARAMS.dna_conc,
                temp_c=THERMO_PARAMS.temp_c,
                max_loop=THERMO_PARAMS.max_loop
            )
            # dG < -9 kcal/mol 或 Tm > 40°C 认为有风险
            has_dimer = result.structure_found and (result.dg < -9000 or result.tm > 40)
            return has_dimer, result.dg / 1000 if result.structure_found else None, result.tm if result.structure_found else None
        except Exception:
            pass
    
    # 回退到简单检测
    return _check_self_dimer_fallback(seq), None, None


# 保持兼容性的别名
def check_self_dimer(seq: str, min_complementary: int = 4) -> bool:
    """检测自身二聚体 (兼容旧接口)"""
    has_dimer, _, _ = check_homodimer(seq)
    return has_dimer


def _check_self_dimer_fallback(seq: str, min_complementary: int = 4) -> bool:
    """内置自身二聚体检测(回退方案)"""
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


def check_heterodimer(seq1: str, seq2: str) -> Tuple[bool, Optional[float], Optional[float]]:
    """
    使用Primer3检测引物间二聚体(异源二聚体)
    
    返回: (是否有二聚体风险, dG值, Tm值)
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    seq1 = re.sub(r'[^ATGC]', '', seq1)
    seq2 = re.sub(r'[^ATGC]', '', seq2)
    
    if len(seq1) < 10 or len(seq2) < 10:
        return False, None, None
    
    if PRIMER3_AVAILABLE:
        try:
            result = primer3.calc_heterodimer(
                seq1, seq2,
                mv_conc=THERMO_PARAMS.mv_conc,
                dv_conc=THERMO_PARAMS.dv_conc,
                dntp_conc=THERMO_PARAMS.dntp_conc,
                dna_conc=THERMO_PARAMS.dna_conc,
                temp_c=THERMO_PARAMS.temp_c,
                max_loop=THERMO_PARAMS.max_loop
            )
            # dG < -9 kcal/mol 或 Tm > 40°C 认为有风险
            has_dimer = result.structure_found and (result.dg < -9000 or result.tm > 40)
            return has_dimer, result.dg / 1000 if result.structure_found else None, result.tm if result.structure_found else None
        except Exception:
            pass
    
    # 回退到简单检测
    return _check_primer_dimer_fallback(seq1, seq2), None, None


# 保持兼容性的别名
def check_primer_dimer(seq1: str, seq2: str, min_complementary: int = 4) -> bool:
    """检测引物二聚体 (兼容旧接口)"""
    has_dimer, _, _ = check_heterodimer(seq1, seq2)
    return has_dimer


def _check_primer_dimer_fallback(seq1: str, seq2: str, min_complementary: int = 4) -> bool:
    """内置引物二聚体检测(回退方案)"""
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


def check_end_stability(seq1: str, seq2: str) -> Tuple[bool, Optional[float]]:
    """
    使用Primer3检测3'端稳定性
    
    返回: (是否稳定, dG值)
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    if PRIMER3_AVAILABLE:
        try:
            result = primer3.calc_end_stability(
                seq1, seq2,
                mv_conc=THERMO_PARAMS.mv_conc,
                dv_conc=THERMO_PARAMS.dv_conc,
                dntp_conc=THERMO_PARAMS.dntp_conc,
                dna_conc=THERMO_PARAMS.dna_conc
            )
            # dG值越负，稳定性越好；但过于稳定可能导致非特异性
            is_good = -10000 < result.dg < -3000  # -10 to -3 kcal/mol
            return is_good, result.dg / 1000
        except Exception:
            pass
    
    return True, None


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


# ==================== KASP人工错配模块 (LGC标准) ====================

"""
KASP人工错配(Deliberate Mismatch)说明:

在KASP检测中，为了增强等位基因特异性(Allele Specificity)，需要在ASP引物的
3'端倒数第2或第3位(n-2或n-3位置)引入人工错配碱基。

原理：
- SNP位点本身只产生一个错配，特异性可能不足
- 引入额外错配可以"削弱"引物与非目标等位基因的结合
- 需要根据SNP碱基的互补强度选择适当的错配

LGC标准错配规则：
- 强互补SNP (G/C)：需要强destabilizing错配
- 弱互补SNP (A/T)：需要中等destabilizing错配
"""

# LGC KASP标准错配规则表
# 格式: {(SNP碱基, 原始n-3碱基): 推荐错配碱基}
# 原理: 根据SNP强度和n-3位置碱基选择最佳destabilizing mismatch

KASP_MISMATCH_RULES = {
    # ==== 强SNP (G或C) - 需要强错配 ====
    # SNP=G时的错配规则
    ('G', 'A'): 'C',  # A→C 强错配
    ('G', 'T'): 'C',  # T→C 强错配
    ('G', 'G'): 'A',  # G→A 强错配
    ('G', 'C'): 'A',  # C→A 强错配
    
    # SNP=C时的错配规则
    ('C', 'A'): 'C',  # A→C 强错配
    ('C', 'T'): 'C',  # T→C 强错配
    ('C', 'G'): 'A',  # G→A 强错配
    ('C', 'C'): 'A',  # C→A 强错配
    
    # ==== 弱SNP (A或T) - 需要中等错配 ====
    # SNP=A时的错配规则
    ('A', 'A'): 'G',  # A→G 中等错配
    ('A', 'T'): 'G',  # T→G 中等错配
    ('A', 'G'): 'T',  # G→T 中等错配
    ('A', 'C'): 'T',  # C→T 中等错配
    
    # SNP=T时的错配规则
    ('T', 'A'): 'G',  # A→G 中等错配
    ('T', 'T'): 'G',  # T→G 中等错配
    ('T', 'G'): 'T',  # G→T 中等错配
    ('T', 'C'): 'T',  # C→T 中等错配
}

# 备用错配规则（当标准规则不适用时）
KASP_FALLBACK_MISMATCH = {
    'A': 'C',  # A的最强destabilizer
    'T': 'C',  # T的最强destabilizer  
    'G': 'A',  # G的最强destabilizer
    'C': 'A',  # C的最强destabilizer
}


def get_kasp_deliberate_mismatch(snp_base: str, original_n3_base: str) -> str:
    """
    获取KASP人工错配碱基 (LGC标准)
    
    参数:
        snp_base: SNP位点的碱基 (引物3'端最后一个碱基)
        original_n3_base: 引物倒数第3位的原始碱基
    
    返回:
        推荐的错配碱基
    
    原理:
        - 强SNP(G/C): 与模板形成3个氢键，需要强错配来destabilize
        - 弱SNP(A/T): 与模板形成2个氢键，需要中等错配
    """
    snp_base = snp_base.upper()
    original_n3_base = original_n3_base.upper()
    
    # 查找LGC标准规则
    key = (snp_base, original_n3_base)
    if key in KASP_MISMATCH_RULES:
        return KASP_MISMATCH_RULES[key]
    
    # 回退到通用规则
    return KASP_FALLBACK_MISMATCH.get(original_n3_base, 'A')


def apply_deliberate_mismatch(core_seq: str, snp_base: str, mismatch_position: int = -3) -> Tuple[str, str, str]:
    """
    对引物序列应用人工错配
    
    参数:
        core_seq: 不含SNP的核心引物序列
        snp_base: SNP碱基 (将添加到3'端)
        mismatch_position: 错配位置 (相对于3'端, 通常为-3)
    
    返回:
        (带错配的完整引物, 原始碱基, 错配碱基)
    """
    if len(core_seq) < abs(mismatch_position):
        # 序列太短，无法引入错配
        return core_seq + snp_base, '', ''
    
    # 计算n-3位置的实际索引（在加入SNP之前）
    # 最终引物 = core_seq + snp_base
    # 所以n-3位置在core_seq中的索引 = len(core_seq) + mismatch_position
    mismatch_idx = len(core_seq) + mismatch_position + 1  # +1因为SNP还未加入
    
    if mismatch_idx < 0 or mismatch_idx >= len(core_seq):
        return core_seq + snp_base, '', ''
    
    original_base = core_seq[mismatch_idx]
    mismatch_base = get_kasp_deliberate_mismatch(snp_base, original_base)
    
    # 如果计算出的错配碱基与原始碱基相同，尝试其他错配
    if mismatch_base == original_base:
        # 使用备用错配
        alternatives = ['A', 'T', 'G', 'C']
        alternatives.remove(original_base)
        mismatch_base = alternatives[0]  # 选择第一个不同的碱基
    
    # 构建带错配的序列
    modified_core = core_seq[:mismatch_idx] + mismatch_base + core_seq[mismatch_idx + 1:]
    final_primer = modified_core + snp_base
    
    return final_primer, original_base, mismatch_base


def get_strong_mismatch(original_base: str) -> str:
    """
    获取强错配碱基 (兼容旧接口)
    
    注意: 这是简化版本，建议使用get_kasp_deliberate_mismatch获得更精确的错配
    """
    return KASP_FALLBACK_MISMATCH.get(original_base.upper(), 'A')


# ============================================================
# 新增函数：强制n-3错配 (按用户命名需求)
# ============================================================

def apply_n3_mismatch(primer_seq: str, snp_base: str) -> Tuple[str, str, str]:
    """
    强制在引物3'端倒数第3位(n-3)引入错配碱基
    
    这是LGC商业标准的核心特异性增强技术。
    
    参数:
        primer_seq: 不含SNP的引物核心序列 (即上游序列)
        snp_base: SNP位点碱基 (将添加到引物3'端)
    
    返回:
        (带错配的完整引物序列, 原始n-3碱基, 替换后的错配碱基)
    
    示例:
        输入: primer_seq="ATGCTCACCACCACTCT", snp_base="A"
        最终引物 = "ATGCTCACCACCACTCT" + "A" + "C" = "ATGCTCACCACCACTCTAC" (21bp)
        n-3位置是从3'端数第3位，即 "T" (index -3 of final primer)
        但n-3位置在core_seq中是最后一位，需要在此处引入错配
        
        实际: n-3 = core_seq[-2]，因为最终引物 = core + SNP
        所以 core[-2] 对应 final_primer[-3]
    """
    if len(primer_seq) < 2:
        # 序列太短，无法引入错配
        return primer_seq + snp_base, '', ''
    
    # n-3位置计算:
    # 最终引物 = core_seq + snp_base
    # final[-1] = snp_base
    # final[-2] = core[-1]  
    # final[-3] = core[-2]  <-- 这是我们要修改的位置
    mismatch_idx = -2  # 在core_seq中的索引
    
    original_base = primer_seq[mismatch_idx].upper()
    
    # 使用LGC标准规则获取错配碱基
    mismatch_base = get_kasp_deliberate_mismatch(snp_base.upper(), original_base)
    
    # 确保错配碱基与原始碱基不同
    if mismatch_base == original_base:
        # 选择最强的destabilizing碱基
        alternatives = {'A': 'C', 'T': 'C', 'G': 'A', 'C': 'A'}
        mismatch_base = alternatives.get(original_base, 'A')
    
    # 构建带错配的引物序列
    modified_core = primer_seq[:mismatch_idx] + mismatch_base + primer_seq[mismatch_idx + 1:]
    final_primer = modified_core + snp_base
    
    return final_primer, original_base, mismatch_base


# ============================================================
# 新增函数：Tm贪婪延伸算法 (5'端动态延伸)
# ============================================================

def optimize_5prime_extension(upstream: str, snp_base: str, config, 
                               target_tm: float = 60.0,
                               min_len: int = 18,
                               max_len: int = 30,
                               apply_mismatch: bool = True) -> List[Dict]:
    """
    Tm贪婪延伸算法：动态向5'端延伸引物以达到目标Tm
    
    复现LGC商业设计中"为了获得更好Tm而多取一个碱基"的行为。
    
    原理:
        1. 从最小长度开始
        2. 计算当前Tm
        3. 如果Tm < 目标值，向5'端延伸1bp
        4. 重复直到Tm达标或达到最大长度
        5. 返回Tm最接近目标的所有候选
    
    参数:
        upstream: SNP上游序列 (不含SNP位点)
        snp_base: SNP碱基
        config: KASPConfig配置对象
        target_tm: 目标Tm值 (默认60°C，LGC标准范围60-62°C)
        min_len: 最小引物长度 (默认18bp)
        max_len: 最大引物长度 (默认30bp)
        apply_mismatch: 是否应用n-3错配 (默认True)
    
    返回:
        候选引物列表，按Tm与目标的接近程度排序
        每个元素: {
            'sequence': 完整引物序列,
            'core_seq': 核心序列(不含SNP),
            'length': 引物长度,
            'tm': 计算的Tm值,
            'tm_diff_from_target': 与目标Tm的差值,
            'gc_content': GC含量,
            'mismatch_info': 错配信息 (如果apply_mismatch=True)
        }
    """
    candidates = []
    
    # 确保有足够的上游序列
    available_len = len(upstream)
    actual_min = max(min_len - 1, 3)  # core_seq长度 = primer_len - 1 (因为要加SNP)
    actual_max = min(max_len - 1, available_len)
    
    if actual_max < actual_min:
        return []
    
    # 贪婪延伸：从最小长度开始，逐步延伸
    best_candidate = None
    best_tm_diff = float('inf')
    
    for core_len in range(actual_min, actual_max + 1):
        # 从上游序列3'端截取
        core_seq = upstream[-core_len:]
        
        # 应用n-3错配
        if apply_mismatch and len(core_seq) >= 2:
            primer_seq, orig_base, mismatch_base = apply_n3_mismatch(core_seq, snp_base)
            mismatch_info = {
                'position': -3,
                'original': orig_base,
                'replacement': mismatch_base,
                'applied': bool(mismatch_base and mismatch_base != orig_base)
            }
        else:
            primer_seq = core_seq + snp_base
            mismatch_info = {'applied': False}
        
        # 计算Tm (只计算核心序列，不含FAM/HEX标签)
        tm_value = calc_tm(primer_seq)
        gc_content = calc_gc_content(primer_seq)
        
        # 计算与目标Tm的差距
        tm_diff = abs(tm_value - target_tm)
        
        candidate = {
            'sequence': primer_seq,
            'core_seq': core_seq,
            'length': len(primer_seq),
            'tm': tm_value,
            'tm_diff_from_target': tm_diff,
            'gc_content': gc_content,
            'mismatch_info': mismatch_info
        }
        candidates.append(candidate)
        
        # 跟踪最佳候选
        if tm_diff < best_tm_diff:
            best_tm_diff = tm_diff
            best_candidate = candidate
        
        # 贪婪优化：如果Tm已经达到或超过目标，可以停止延伸
        # 但继续搜索几个以找到最优
        if tm_value >= target_tm and len(candidates) >= 3:
            # 再延伸2-3bp看看是否有更好的
            extra_search = 3
            for extra_len in range(core_len + 1, min(core_len + extra_search + 1, actual_max + 1)):
                extra_core = upstream[-extra_len:]
                if apply_mismatch and len(extra_core) >= 2:
                    extra_primer, e_orig, e_mismatch = apply_n3_mismatch(extra_core, snp_base)
                    e_mismatch_info = {
                        'position': -3,
                        'original': e_orig,
                        'replacement': e_mismatch,
                        'applied': bool(e_mismatch and e_mismatch != e_orig)
                    }
                else:
                    extra_primer = extra_core + snp_base
                    e_mismatch_info = {'applied': False}
                
                extra_tm = calc_tm(extra_primer)
                extra_diff = abs(extra_tm - target_tm)
                
                candidates.append({
                    'sequence': extra_primer,
                    'core_seq': extra_core,
                    'length': len(extra_primer),
                    'tm': extra_tm,
                    'tm_diff_from_target': extra_diff,
                    'gc_content': calc_gc_content(extra_primer),
                    'mismatch_info': e_mismatch_info
                })
            break
    
    # 按Tm与目标的接近程度排序
    candidates.sort(key=lambda x: (x['tm_diff_from_target'], -x['length']))
    
    return candidates


def design_optimal_asp_primer(upstream: str, snp_base: str, config,
                               target_tm_range: Tuple[float, float] = (60.0, 62.0)) -> Optional[Dict]:
    """
    设计最优ASP引物：结合n-3错配和Tm贪婪延伸
    
    这是对LGC商业设计的完整复现，包括：
    1. 强制n-3错配增强特异性
    2. 动态5'延伸达到最佳Tm (60-62°C)
    
    参数:
        upstream: SNP上游序列
        snp_base: SNP碱基
        config: KASPConfig配置
        target_tm_range: 目标Tm范围 (默认60-62°C)
    
    返回:
        最优引物信息字典，或None如果无法设计
    """
    target_tm = (target_tm_range[0] + target_tm_range[1]) / 2  # 目标Tm取中间值
    
    # 使用Tm贪婪延伸算法获取候选
    candidates = optimize_5prime_extension(
        upstream=upstream,
        snp_base=snp_base,
        config=config,
        target_tm=target_tm,
        min_len=getattr(config, 'MIN_PRIMER_LEN', 18),
        max_len=getattr(config, 'MAX_PRIMER_LEN', 28),
        apply_mismatch=True
    )
    
    if not candidates:
        return None
    
    # 筛选Tm在目标范围内的候选
    in_range = [c for c in candidates 
                if target_tm_range[0] <= c['tm'] <= target_tm_range[1]]
    
    if in_range:
        # 选择Tm最接近目标中心的
        return min(in_range, key=lambda x: x['tm_diff_from_target'])
    else:
        # 没有完全在范围内的，返回最接近的
        return candidates[0]


def evaluate_primer_quality(seq: str, config=None) -> Dict:
    """综合评估引物质量"""
    if config is None:
        config = KASPConfig()
    
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)
    
    # 使用Primer3进行热力学分析
    tm_value = calc_tm(seq)
    hairpin_result = check_hairpin(seq)
    homodimer_result = check_homodimer(seq)
    
    # 处理返回值 (新API返回元组)
    has_hairpin = hairpin_result[0] if isinstance(hairpin_result, tuple) else hairpin_result
    hairpin_tm = hairpin_result[1] if isinstance(hairpin_result, tuple) and len(hairpin_result) > 1 else None
    
    has_self_dimer = homodimer_result[0] if isinstance(homodimer_result, tuple) else homodimer_result
    homodimer_dg = homodimer_result[1] if isinstance(homodimer_result, tuple) and len(homodimer_result) > 1 else None
    
    result = {
        'sequence': seq,
        'length': len(seq),
        'gc_content': calc_gc_content(seq),
        'tm': tm_value,
        'has_hairpin': has_hairpin,
        'hairpin_tm': hairpin_tm,
        'has_self_dimer': has_self_dimer,
        'homodimer_dg': homodimer_dg,
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


def design_kasp_common_primer_with_primer3(downstream: str, config: KASPConfig,
                                            min_distance: int = 30, max_distance: int = 80,
                                            target_tm: float = None) -> List[Dict]:
    """
    使用Primer3设计KASP通用反向引物 - 优化Tm平衡版本
    
    参数:
        downstream: SNP下游序列
        config: KASP配置
        min_distance: 距SNP最小距离
        max_distance: 距SNP最大距离
        target_tm: 目标Tm值（用于与ASP引物匹配），None则使用配置默认值
    
    返回:
        反向引物候选列表（按Tm匹配度排序）
    """
    if not PRIMER3_AVAILABLE:
        return []
    
    if len(downstream) < max_distance + config.MIN_PRIMER_LEN:
        max_distance = len(downstream) - config.MIN_PRIMER_LEN
    
    if max_distance < min_distance:
        return []
    
    # 如果没有指定目标Tm，使用Common引物目标Tm (60-62°C)
    if target_tm is None:
        target_tm = getattr(config, 'COMMON_TARGET_TM', 61.0)
    
    # Common引物可以使用更长的长度来达到目标Tm
    common_max_len = getattr(config, 'COMMON_MAX_LEN', 30)
    
    candidates = []
    
    # 使用Primer3的设计引擎只设计右引物
    seq_args = {
        'SEQUENCE_ID': 'kasp_common',
        'SEQUENCE_TEMPLATE': downstream,
        'SEQUENCE_FORCE_LEFT_START': 0,
    }
    
    global_args = {
        'PRIMER_TASK': 'pick_primer_list',
        'PRIMER_PICK_LEFT_PRIMER': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_NUM_RETURN': 30,  # 返回更多候选以便筛选
        
        'PRIMER_MIN_SIZE': config.MIN_PRIMER_LEN,
        'PRIMER_OPT_SIZE': 22,  # Common引物最优长度稍长
        'PRIMER_MAX_SIZE': common_max_len,  # 允许更长以达到Tm目标
        
        # Tm参数 - 针对Common引物优化
        'PRIMER_MIN_TM': target_tm - 3.0,
        'PRIMER_OPT_TM': target_tm,
        'PRIMER_MAX_TM': target_tm + 3.0,
        
        'PRIMER_MIN_GC': config.MIN_GC,
        'PRIMER_MAX_GC': config.MAX_GC,
        
        'PRIMER_SALT_MONOVALENT': THERMO_PARAMS.kasp_mv_conc,
        'PRIMER_SALT_DIVALENT': THERMO_PARAMS.kasp_dv_conc,
        'PRIMER_DNTP_CONC': THERMO_PARAMS.kasp_dntp_conc,
        'PRIMER_DNA_CONC': THERMO_PARAMS.kasp_dna_conc,
        
        'PRIMER_MAX_SELF_ANY': 8,
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_MAX_HAIRPIN_TH': 47.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_GC_CLAMP': 1,
        
        'PRIMER_PRODUCT_SIZE_RANGE': [[min_distance + config.MIN_PRIMER_LEN, 
                                       max_distance + common_max_len]],
    }
    
    try:
        results = primer3.design_primers(seq_args, global_args)
        num_returned = results.get('PRIMER_RIGHT_NUM_RETURNED', 0)
        
        for i in range(num_returned):
            try:
                seq = results.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
                pos = results.get(f'PRIMER_RIGHT_{i}', [0, 0])
                tm = results.get(f'PRIMER_RIGHT_{i}_TM', 0)
                gc = results.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0)
                
                if not seq:
                    continue
                
                # 计算距离SNP的距离
                rev_distance = pos[0] - pos[1] + 1
                
                # 检查距离是否在范围内
                if rev_distance < min_distance or rev_distance > max_distance:
                    continue
                
                # 检测二级结构
                has_hairpin, hairpin_tm = check_hairpin(seq)
                has_dimer, dimer_dg, _ = check_homodimer(seq)
                
                # 计算与目标Tm的偏差（用于排序）
                tm_deviation = abs(tm - target_tm)
                
                candidate = {
                    'sequence': seq,
                    'position': pos,
                    'distance': rev_distance,
                    'tm': round(tm, 1),
                    'gc': round(gc, 1),
                    'has_hairpin': has_hairpin,
                    'hairpin_tm': hairpin_tm,
                    'has_dimer': has_dimer,
                    'dimer_dg': dimer_dg,
                    'penalty': results.get(f'PRIMER_RIGHT_{i}_PENALTY', 0),
                    'tm_deviation': tm_deviation,
                    'target_tm': target_tm
                }
                candidates.append(candidate)
                
            except Exception:
                continue
        
        # 按Tm偏差排序（优先选择Tm接近目标的引物）
        candidates.sort(key=lambda x: (x['tm_deviation'], x['penalty']))
        
    except Exception:
        pass
    
    return candidates


# ==================== 智能失败诊断模块 (Smart Failure Feedback) ====================

@dataclass
class DesignFailureDiagnosis:
    """引物设计失败诊断结果"""
    success: bool                      # 是否成功设计
    failure_reason: str                # 主要失败原因
    suggestions: List[str]             # 优化建议列表
    sequence_analysis: Dict            # 序列分析详情
    severity: str                      # 问题严重程度: 'low', 'medium', 'high', 'critical'


def diagnose_design_failure(upstream: str, downstream: str, allele1: str, allele2: str,
                            config: KASPConfig = None) -> DesignFailureDiagnosis:
    """
    智能失败诊断 - 分析为什么引物设计失败并给出具体建议
    
    参数:
        upstream: SNP上游序列
        downstream: SNP下游序列
        allele1, allele2: 两个等位基因
        config: KASP配置
    
    返回:
        DesignFailureDiagnosis 诊断结果对象
    """
    if config is None:
        config = KASPConfig()
    
    suggestions = []
    issues_found = []
    severity = 'low'
    
    # === 1. 分析目标区域 (SNP上下游50bp) ===
    target_upstream = upstream[-50:] if len(upstream) >= 50 else upstream
    target_downstream = downstream[:50] if len(downstream) >= 50 else downstream
    target_region = target_upstream + downstream[:80]  # 引物设计核心区域
    
    # === 2. GC含量分析 ===
    target_gc = calc_gc_content(target_region)
    upstream_gc = calc_gc_content(target_upstream)
    downstream_gc = calc_gc_content(target_downstream)
    
    gc_issues = []
    
    # 极低GC含量检测 (< 25%) - 关键问题
    if target_gc < 25:
        gc_issues.append(f"目标区域GC含量极低 ({target_gc:.1f}% < 25%)")
        issues_found.append("gc_extremely_low")
        severity = 'critical'
        suggestions.append("🔴 该位点GC含量极低（<25%），常规引物难以结合，强烈建议重新选择附近的SNP位点")
        suggestions.append("💡 尝试在该SNP上下游100-200bp范围内寻找GC含量更高的替代SNP")
    elif target_gc < 30:
        gc_issues.append(f"目标区域GC含量偏低 ({target_gc:.1f}%)")
        issues_found.append("gc_low")
        if severity != 'critical':
            severity = 'high'
        suggestions.append("🟠 GC含量偏低，可能导致引物结合力不足")
        suggestions.append("💡 可尝试增加引物长度(28-32bp)以提高Tm值")
    
    # 上游/下游GC不平衡
    gc_diff = abs(upstream_gc - downstream_gc)
    if gc_diff > 25:
        gc_issues.append(f"上下游GC含量差异大 (上游{upstream_gc:.1f}% vs 下游{downstream_gc:.1f}%)")
        issues_found.append("gc_imbalance")
        suggestions.append("⚠️ 上下游GC含量差异较大，可能导致ASP和Common引物Tm难以平衡")
    
    # === 3. 序列复杂度分析 ===
    complexity = analyze_sequence_complexity(target_region)
    
    if complexity['low_complexity']:
        issues_found.append("low_complexity")
        if severity not in ['critical', 'high']:
            severity = 'medium'
        suggestions.append("⚠️ 序列含有低复杂度区域（连续重复碱基），可能导致非特异性扩增")
    
    if complexity['complexity_score'] < 50:
        issues_found.append("simple_sequence")
        suggestions.append("⚠️ 序列复杂度较低，引物可能在基因组中有多个结合位点")
    
    # === 4. 发夹结构风险分析 ===
    # 检查上游序列（ASP引物区域）
    asp_region = upstream[-30:] if len(upstream) >= 30 else upstream
    has_hairpin_asp, hairpin_tm_asp = check_hairpin(asp_region)
    
    if has_hairpin_asp:
        issues_found.append("hairpin_risk")
        if severity not in ['critical']:
            severity = 'high'
        suggestions.append(f"🔺 ASP引物区域存在发夹结构风险 (Tm: {hairpin_tm_asp}°C)")
        suggestions.append("💡 可尝试向5'端延伸引物，避开发夹结构区域")
    
    # === 5. 重复序列检测 ===
    has_repeat, repeat_issues = check_wheat_repeat_sequences(target_region)
    if has_repeat:
        issues_found.append("repeat_sequences")
        if severity not in ['critical', 'high']:
            severity = 'medium'
        for issue in repeat_issues:
            suggestions.append(f"🔁 检测到重复序列特征: {issue}")
        suggestions.append("💡 该区域可能位于转座子或SSR区域，引物特异性可能受影响")
    
    # === 6. 序列长度检查 ===
    if len(upstream) < 25:
        issues_found.append("upstream_too_short")
        severity = 'critical'
        suggestions.append("🔴 上游序列过短（<25bp），无法设计有效的ASP引物")
        suggestions.append("💡 请提供更长的上游侧翼序列（建议>50bp）")
    
    if len(downstream) < 60:
        issues_found.append("downstream_too_short")
        if severity != 'critical':
            severity = 'high'
        suggestions.append("🟠 下游序列较短，Common引物设计空间受限")
        suggestions.append("💡 请提供更长的下游侧翼序列（建议>100bp）")
    
    # === 7. SNP类型分析 ===
    snp_pair = frozenset([allele1.upper(), allele2.upper()])
    
    # 转换/颠换分析
    transitions = [frozenset(['A', 'G']), frozenset(['C', 'T'])]
    transversions = [frozenset(['A', 'C']), frozenset(['A', 'T']), 
                     frozenset(['G', 'C']), frozenset(['G', 'T'])]
    
    if snp_pair in transitions:
        snp_type = "转换(Transition)"
    elif snp_pair in transversions:
        snp_type = "颠换(Transversion)"
    else:
        snp_type = "未知"
    
    # AT/AT SNP特别难设计
    if snp_pair == frozenset(['A', 'T']):
        issues_found.append("at_snp")
        suggestions.append("⚠️ A/T SNP的两个等位基因引物Tm差异可能较大")
    
    # === 8. 综合诊断 ===
    sequence_analysis = {
        'target_gc': target_gc,
        'upstream_gc': upstream_gc,
        'downstream_gc': downstream_gc,
        'gc_difference': gc_diff,
        'complexity_score': complexity['complexity_score'],
        'has_low_complexity': complexity['low_complexity'],
        'has_hairpin': has_hairpin_asp,
        'hairpin_tm': hairpin_tm_asp,
        'has_repeat': has_repeat,
        'upstream_length': len(upstream),
        'downstream_length': len(downstream),
        'snp_type': snp_type,
        'snp_bases': f"[{allele1}/{allele2}]",
        'issues_found': issues_found
    }
    
    # 生成主要失败原因
    if 'gc_extremely_low' in issues_found:
        failure_reason = f"目标区域GC含量极低 ({target_gc:.1f}%)，无法设计出满足Tm要求的引物"
    elif 'upstream_too_short' in issues_found:
        failure_reason = "上游序列过短，ASP引物设计空间不足"
    elif 'gc_low' in issues_found and 'hairpin_risk' in issues_found:
        failure_reason = "GC含量偏低且存在发夹结构风险，双重限制导致设计失败"
    elif 'low_complexity' in issues_found or 'repeat_sequences' in issues_found:
        failure_reason = "该区域序列复杂度低或存在重复序列，引物设计空间受限"
    elif 'hairpin_risk' in issues_found:
        failure_reason = "引物结合区域存在强发夹结构，无法找到合适的引物"
    elif 'gc_low' in issues_found:
        failure_reason = f"GC含量偏低 ({target_gc:.1f}%)，难以达到目标Tm值"
    else:
        failure_reason = "综合因素导致无法找到满足所有条件的引物组合"
        suggestions.append("💡 尝试放宽Tm范围（如52-65°C）或增加最大引物长度（如35bp）")
    
    # 添加通用建议
    if not suggestions:
        suggestions = [
            "💡 检查输入序列是否正确",
            "💡 尝试提供更长的侧翼序列",
            "💡 考虑选择附近的替代SNP位点"
        ]
    
    return DesignFailureDiagnosis(
        success=False,
        failure_reason=failure_reason,
        suggestions=suggestions,
        sequence_analysis=sequence_analysis,
        severity=severity
    )


def format_diagnosis_for_display(diagnosis: DesignFailureDiagnosis) -> str:
    """将诊断结果格式化为用户友好的显示文本"""
    severity_icons = {
        'critical': '🔴',
        'high': '🟠',
        'medium': '🟡',
        'low': '🟢'
    }
    
    severity_labels = {
        'critical': '严重',
        'high': '较高',
        'medium': '中等',
        'low': '轻微'
    }
    
    icon = severity_icons.get(diagnosis.severity, '⚪')
    label = severity_labels.get(diagnosis.severity, '未知')
    
    output = f"""
### {icon} 设计失败诊断报告

**问题严重程度**: {label}

**主要原因**: {diagnosis.failure_reason}

---

#### 📊 序列分析结果

| 指标 | 数值 | 状态 |
|------|------|------|
| 目标区域GC% | {diagnosis.sequence_analysis['target_gc']:.1f}% | {'⚠️ 偏低' if diagnosis.sequence_analysis['target_gc'] < 35 else '✓ 正常'} |
| 上游序列GC% | {diagnosis.sequence_analysis['upstream_gc']:.1f}% | {'⚠️' if diagnosis.sequence_analysis['upstream_gc'] < 30 else '✓'} |
| 下游序列GC% | {diagnosis.sequence_analysis['downstream_gc']:.1f}% | {'⚠️' if diagnosis.sequence_analysis['downstream_gc'] < 30 else '✓'} |
| 序列复杂度 | {diagnosis.sequence_analysis['complexity_score']:.0f}/100 | {'⚠️ 偏低' if diagnosis.sequence_analysis['complexity_score'] < 60 else '✓ 正常'} |
| 发夹结构风险 | {'是' if diagnosis.sequence_analysis['has_hairpin'] else '否'} | {'⚠️' if diagnosis.sequence_analysis['has_hairpin'] else '✓'} |
| 重复序列 | {'检测到' if diagnosis.sequence_analysis['has_repeat'] else '未检测到'} | {'⚠️' if diagnosis.sequence_analysis['has_repeat'] else '✓'} |
| 上游序列长度 | {diagnosis.sequence_analysis['upstream_length']}bp | {'⚠️ 过短' if diagnosis.sequence_analysis['upstream_length'] < 30 else '✓'} |
| 下游序列长度 | {diagnosis.sequence_analysis['downstream_length']}bp | {'⚠️ 过短' if diagnosis.sequence_analysis['downstream_length'] < 80 else '✓'} |
| SNP类型 | {diagnosis.sequence_analysis['snp_type']} {diagnosis.sequence_analysis['snp_bases']} | - |

---

#### 💡 优化建议

"""
    for i, suggestion in enumerate(diagnosis.suggestions, 1):
        output += f"{i}. {suggestion}\n"
    
    return output


def design_kasp_primers_multi(upstream: str, downstream: str, allele1: str, allele2: str, 
                              config: KASPConfig = None, num_schemes: int = 5,
                              _is_rescue_mode: bool = False) -> List[Dict]:
    """
    设计多套KASP引物方案 - 优化版 (支持Primer3、人工错配、救援模式、Tm平衡)
    
    核心优化:
    1. LGC标准人工错配(Deliberate Mismatch) - 增强等位基因特异性
    2. 救援模式(Rescue Mode) - 处理AT-rich序列
    3. Tm平衡优化 - ASP与Common引物Tm匹配
    
    参数:
        upstream: SNP上游序列
        downstream: SNP下游序列  
        allele1, allele2: 两个等位基因碱基
        config: KASP配置
        num_schemes: 需要返回的方案数
        _is_rescue_mode: 内部参数，标记是否处于救援模式
    
    返回:
        引物方案列表，按质量评分排序
    """
    if config is None:
        config = KASPConfig()
    
    # 检查序列长度是否足够
    if len(upstream) < config.MIN_PRIMER_LEN:
        return []  # 上游序列太短
    
    if len(downstream) < config.REV_MIN_DISTANCE + config.MIN_PRIMER_LEN:
        return []  # 下游序列太短
    
    all_schemes = []
    seen_signatures = set()  # 在线去重：记录已添加的方案签名
    
    # === 第一步：设计ASP引物候选，计算平均Tm用于Common引物匹配 ===
    asp_tm_values = []  # 收集ASP引物Tm值用于计算目标Tm
    
    # === Tm贪婪延伸优化 (LGC商业设计标准) ===
    # 如果启用，使用动态延伸算法找到最佳Tm的引物长度
    asp_target_tm_min = getattr(config, 'ASP_TARGET_TM_MIN', 60.0)
    asp_target_tm_max = getattr(config, 'ASP_TARGET_TM_MAX', 62.0)
    use_tm_greedy = getattr(config, 'TM_GREEDY_EXTENSION', True)
    
    # 收集所有ASP候选（含Tm信息）用于智能筛选
    asp_candidates = []
    
    # 生成不同长度的ASP引物
    for primer_len in range(config.MIN_PRIMER_LEN, min(config.MAX_PRIMER_LEN + 1, len(upstream) + 1)):
        core_seq = upstream[-(primer_len - 1):]
        
        if len(core_seq) < config.MIN_PRIMER_LEN - 1:
            continue
        
        # 小麦模式：预筛选GC含量
        if config.WHEAT_MODE and config.WHEAT_STRICT_GC:
            core_gc = calc_gc_content(core_seq)
            if core_gc > 65 or core_gc < 30:
                continue
        
        # === 使用LGC标准人工错配 ===
        for mismatch_pos in config.MISMATCH_POSITIONS:
            if abs(mismatch_pos) >= len(core_seq):
                continue
            
            # 应用人工错配到两个等位基因引物
            # 注意：mismatch_pos是相对于最终引物3'端的位置
            fwd_allele1, orig_base1, mismatch_base1 = apply_deliberate_mismatch(
                core_seq, allele1, mismatch_position=mismatch_pos
            )
            fwd_allele2, orig_base2, mismatch_base2 = apply_deliberate_mismatch(
                core_seq, allele2, mismatch_position=mismatch_pos
            )
            
            # 验证错配是否有效（两个引物应该有相同的错配位置和变化）
            if not mismatch_base1 or mismatch_base1 == orig_base1:
                continue
            
            # 添加FAM/HEX尾巴
            fwd_with_fam = config.FAM_TAIL + fwd_allele1
            fwd_with_hex = config.HEX_TAIL + fwd_allele2
            
            # === ASP引物评估 (只计算核心序列的Tm，不含标签) ===
            eval1 = evaluate_primer_quality(fwd_allele1, config)
            eval2 = evaluate_primer_quality(fwd_allele2, config)
            
            # 质量过低的直接跳过
            if eval1['score'] < 40 or eval2['score'] < 40:
                continue
            
            # ASP引物间的Tm差异
            asp_tm_diff = abs(eval1['tm'] - eval2['tm'])
            
            # Tm差异过大直接跳过
            if asp_tm_diff > config.MAX_TM_DIFF + 1:
                continue
            
            # 收集ASP的Tm用于Common引物设计
            asp_avg_tm = (eval1['tm'] + eval2['tm']) / 2
            asp_tm_values.append(asp_avg_tm)
            
            # === Tm贪婪延伸检查 ===
            # 如果Tm低于目标范围下限，尝试向5'端延伸
            tm_optimized = False
            if use_tm_greedy and asp_avg_tm < asp_target_tm_min:
                # 计算还需要多少Tm提升
                tm_deficit = asp_target_tm_min - asp_avg_tm
                
                # 尝试延伸以达到目标Tm
                extended_len = primer_len
                max_extend = min(config.MAX_PRIMER_LEN + 5, len(upstream))  # 允许额外延伸5bp
                
                while extended_len < max_extend and asp_avg_tm < asp_target_tm_min:
                    extended_len += 1
                    ext_core_seq = upstream[-(extended_len - 1):]
                    
                    # 重新应用错配
                    ext_fwd1, ext_orig1, ext_mis1 = apply_deliberate_mismatch(
                        ext_core_seq, allele1, mismatch_position=mismatch_pos
                    )
                    ext_fwd2, ext_orig2, ext_mis2 = apply_deliberate_mismatch(
                        ext_core_seq, allele2, mismatch_position=mismatch_pos
                    )
                    
                    if not ext_mis1 or ext_mis1 == ext_orig1:
                        continue
                    
                    # 重新计算Tm
                    ext_eval1 = evaluate_primer_quality(ext_fwd1, config)
                    ext_eval2 = evaluate_primer_quality(ext_fwd2, config)
                    ext_avg_tm = (ext_eval1['tm'] + ext_eval2['tm']) / 2
                    
                    # 如果Tm达到目标，使用延伸后的引物
                    if ext_avg_tm >= asp_target_tm_min - 0.5:
                        fwd_allele1, fwd_allele2 = ext_fwd1, ext_fwd2
                        eval1, eval2 = ext_eval1, ext_eval2
                        orig_base1, mismatch_base1 = ext_orig1, ext_mis1
                        core_seq = ext_core_seq
                        asp_avg_tm = ext_avg_tm
                        asp_tm_diff = abs(eval1['tm'] - eval2['tm'])
                        fwd_with_fam = config.FAM_TAIL + fwd_allele1
                        fwd_with_hex = config.HEX_TAIL + fwd_allele2
                        tm_optimized = True
                        break
            
            # 用于错配信息记录
            original_base = orig_base1
            mismatch_base = mismatch_base1
            
            # === 设计Common引物，目标Tm与ASP匹配 ===
            # Common引物目标Tm略高于ASP，因为ASP有竞争反应
            common_target_tm = asp_avg_tm + 0.5  # 略高0.5°C
            common_target_tm = max(
                getattr(config, 'COMMON_TARGET_TM', 61.0) - 2,
                min(getattr(config, 'COMMON_TARGET_TM', 61.0) + 2, common_target_tm)
            )
            
            # 使用Primer3设计Common引物（带Tm目标）
            primer3_common_candidates = []
            if PRIMER3_AVAILABLE:
                primer3_common_candidates = design_kasp_common_primer_with_primer3(
                    downstream, config,
                    min_distance=config.REV_MIN_DISTANCE,
                    max_distance=config.REV_MAX_DISTANCE,
                    target_tm=common_target_tm
                )
            
            # === 匹配Common引物 ===
            if primer3_common_candidates:
                for common_cand in primer3_common_candidates[:10]:
                    rev_seq = common_cand['sequence']
                    rev_dist = common_cand['distance']
                    
                    eval_rev = evaluate_primer_quality(rev_seq, config)
                    
                    if eval_rev['score'] < 40:
                        continue
                    
                    # 检查引物二聚体
                    has_dimer_1, dimer_dg_1, _ = check_heterodimer(fwd_allele1, rev_seq)
                    has_dimer_2, dimer_dg_2, _ = check_heterodimer(fwd_allele2, rev_seq)
                    has_dimer = has_dimer_1 or has_dimer_2
                    
                    product_size = len(upstream) + 1 + rev_dist + len(rev_seq)
                    
                    # === Tm平衡评分 ===
                    # ASP与Common的Tm差异
                    asp_common_tm_diff = abs(asp_avg_tm - eval_rev['tm'])
                    asp_common_max_diff = getattr(config, 'ASP_COMMON_TM_DIFF_MAX', 3.0)
                    
                    # 综合评分
                    avg_fwd_score = (eval1['score'] + eval2['score']) / 2
                    total_score = (avg_fwd_score * 0.35 + eval_rev['score'] * 0.35)
                    
                    # ASP间Tm匹配评分
                    if asp_tm_diff <= 0.5:
                        total_score += 12
                    elif asp_tm_diff <= 1.0:
                        total_score += 8
                    elif asp_tm_diff <= 2.0:
                        total_score += 4
                    else:
                        total_score -= 8
                    
                    # ASP-Common Tm匹配评分（新增）
                    if asp_common_tm_diff <= 1.0:
                        total_score += 10
                    elif asp_common_tm_diff <= 2.0:
                        total_score += 5
                    elif asp_common_tm_diff > asp_common_max_diff:
                        total_score -= 10
                    
                    if has_dimer:
                        total_score -= 15
                    
                    # Primer3优化加分
                    total_score += 5
                    
                    # === Tm贪婪延伸优化加分 ===
                    # 如果Tm在目标范围内，额外加分
                    if asp_target_tm_min <= asp_avg_tm <= asp_target_tm_max:
                        total_score += 8  # Tm达到LGC标准范围
                    
                    wheat_issues = []
                    wheat_details = {}
                    
                    if config.WHEAT_MODE:
                        wheat_bonus, wheat_issues, wheat_details = evaluate_kasp_wheat_specificity(
                            upstream, downstream, fwd_allele1, rev_seq, config
                        )
                        total_score += wheat_bonus
                        amplicon_status, amplicon_bonus = check_amplicon_length_kasp(product_size)
                        total_score += amplicon_bonus
                        wheat_details['amplicon_status'] = amplicon_status
                    
                    total_score = max(0, min(100, total_score))
                    
                    is_usable = True
                    if config.WHEAT_MODE:
                        is_usable = (
                            total_score >= 50 and
                            product_size <= 120 and
                            30 <= eval1['gc_content'] <= 65 and
                            30 <= eval_rev['gc_content'] <= 65 and
                            not has_dimer and
                            asp_tm_diff <= config.MAX_TM_DIFF
                        )
                    else:
                        is_usable = (
                            total_score >= 45 and
                            not has_dimer and
                            asp_tm_diff <= config.MAX_TM_DIFF
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
                        'deliberate_mismatch_info': {
                            'snp_bases': (allele1, allele2),
                            'position': mismatch_pos,
                            'original': original_base,
                            'replacement': mismatch_base,
                            'rule': 'LGC_STANDARD'
                        },
                        'eval_fwd1': eval1,
                        'eval_fwd2': eval2,
                        'eval_rev': eval_rev,
                        'tm_diff': asp_tm_diff,
                        'asp_common_tm_diff': asp_common_tm_diff,
                        'has_dimer': has_dimer,
                        'heterodimer_dg': min(dimer_dg_1 or 0, dimer_dg_2 or 0) if (dimer_dg_1 or dimer_dg_2) else None,
                        'product_size': product_size,
                        'rev_distance': rev_dist,
                        'total_score': total_score,
                        'is_usable': is_usable,
                        'wheat_mode': config.WHEAT_MODE,
                        'wheat_issues': wheat_issues,
                        'wheat_details': wheat_details,
                        'primer3_designed': True,
                        'rescue_mode': _is_rescue_mode,
                        'tm_greedy_optimized': tm_optimized,  # 标记是否经过Tm贪婪延伸优化
                        'asp_target_tm_range': (asp_target_tm_min, asp_target_tm_max)
                    }
                    # 在线去重：检查签名是否已存在
                    signature = (fwd_allele1, fwd_allele2, rev_seq, mismatch_pos)
                    if signature not in seen_signatures:
                        seen_signatures.add(signature)
                        all_schemes.append(scheme)
            
            # === 回退：手动搜索反向引物 ===
            common_max_len = getattr(config, 'COMMON_MAX_LEN', config.MAX_PRIMER_LEN)
            max_rev_dist = min(config.REV_MAX_DISTANCE + 1, len(downstream) - config.MIN_PRIMER_LEN + 1)
            
            if max_rev_dist <= config.REV_MIN_DISTANCE:
                continue
                
            for rev_dist in range(config.REV_MIN_DISTANCE, max_rev_dist):
                for rev_len in range(config.MIN_PRIMER_LEN, min(common_max_len + 1, len(downstream) - rev_dist + 1)):
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
                    
                    # Tm平衡检查
                    asp_common_tm_diff = abs(asp_avg_tm - eval_rev['tm'])
                    asp_common_max_diff = getattr(config, 'ASP_COMMON_TM_DIFF_MAX', 3.0)
                    
                    # 检查引物二聚体
                    has_dimer = (check_primer_dimer(fwd_allele1, rev_seq) or 
                                check_primer_dimer(fwd_allele2, rev_seq))
                    
                    # 计算产物大小
                    product_size = len(upstream) + 1 + rev_dist + rev_len
                    
                    # === 综合评分 ===
                    avg_fwd_score = (eval1['score'] + eval2['score']) / 2
                    total_score = (avg_fwd_score * 0.35 + eval_rev['score'] * 0.35)
                    
                    # ASP间Tm匹配评分
                    if asp_tm_diff <= 0.5:
                        total_score += 12
                    elif asp_tm_diff <= 1.0:
                        total_score += 8
                    elif asp_tm_diff <= 2.0:
                        total_score += 4
                    else:
                        total_score -= 8
                    
                    # ASP-Common Tm匹配评分
                    if asp_common_tm_diff <= 1.0:
                        total_score += 10
                    elif asp_common_tm_diff <= 2.0:
                        total_score += 5
                    elif asp_common_tm_diff > asp_common_max_diff:
                        total_score -= 10
                    
                    # 二聚体惩罚
                    if has_dimer:
                        total_score -= 15
                    
                    # 小麦特异性评分
                    wheat_issues = []
                    wheat_details = {}
                    
                    if config.WHEAT_MODE:
                        wheat_bonus, wheat_issues, wheat_details = evaluate_kasp_wheat_specificity(
                            upstream, downstream, fwd_allele1, rev_seq, config
                        )
                        total_score += wheat_bonus
                        
                        amplicon_status, amplicon_bonus = check_amplicon_length_kasp(product_size)
                        total_score += amplicon_bonus
                        wheat_details['amplicon_status'] = amplicon_status
                    else:
                        if config.PRODUCT_MIN <= product_size <= config.PRODUCT_MAX:
                            total_score += 5
                    
                    total_score = max(0, min(100, total_score))
                    
                    # 判断是否可用
                    is_usable = True
                    if config.WHEAT_MODE:
                        is_usable = (
                            total_score >= 50 and
                            product_size <= 120 and
                            30 <= eval1['gc_content'] <= 65 and
                            30 <= eval_rev['gc_content'] <= 65 and
                            not has_dimer and
                            asp_tm_diff <= config.MAX_TM_DIFF
                        )
                    else:
                        is_usable = (
                            total_score >= 45 and
                            not has_dimer and
                            asp_tm_diff <= config.MAX_TM_DIFF
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
                        'deliberate_mismatch_info': {
                            'snp_bases': (allele1, allele2),
                            'position': mismatch_pos,
                            'original': original_base,
                            'replacement': mismatch_base,
                            'rule': 'LGC_STANDARD'
                        },
                        'eval_fwd1': eval1,
                        'eval_fwd2': eval2,
                        'eval_rev': eval_rev,
                        'tm_diff': asp_tm_diff,
                        'asp_common_tm_diff': asp_common_tm_diff,
                        'has_dimer': has_dimer,
                        'product_size': product_size,
                        'rev_distance': rev_dist,
                        'total_score': total_score,
                        'is_usable': is_usable,
                        'wheat_mode': config.WHEAT_MODE,
                        'wheat_issues': wheat_issues,
                        'wheat_details': wheat_details,
                        'rescue_mode': _is_rescue_mode,
                        'tm_greedy_optimized': tm_optimized,  # 标记是否经过Tm贪婪延伸优化
                        'asp_target_tm_range': (asp_target_tm_min, asp_target_tm_max)
                    }
                    # 在线去重：检查签名是否已存在
                    signature = (fwd_allele1, fwd_allele2, rev_seq, mismatch_pos)
                    if signature not in seen_signatures:
                        seen_signatures.add(signature)
                        all_schemes.append(scheme)
    
    # === 救援模式：如果标准模式没有找到方案，启用救援模式 ===
    if not all_schemes and config.RESCUE_MODE_ENABLED and not _is_rescue_mode:
        rescue_config = config.get_rescue_config()
        rescue_schemes = design_kasp_primers_multi(
            upstream, downstream, allele1, allele2,
            config=rescue_config,
            num_schemes=num_schemes,
            _is_rescue_mode=True  # 标记为救援模式，防止无限递归
        )
        if rescue_schemes:
            # 标记为救援模式设计的引物
            for scheme in rescue_schemes:
                scheme['rescue_mode'] = True
                scheme['rescue_note'] = "⚠️ 该引物由救援模式设计（放宽参数），建议优先考虑标准模式引物"
            return rescue_schemes
    
    # === 如果救援模式也失败，返回空列表（诊断在外层处理）===
    if not all_schemes:
        return []
    
    # 按评分排序
    if config.WHEAT_MODE:
        all_schemes.sort(key=lambda x: (x.get('is_usable', False), x['total_score']), reverse=True)
    else:
        all_schemes.sort(key=lambda x: x['total_score'], reverse=True)
    
    # === 最终去重（双重保险）：确保返回的方案都是唯一的 ===
    unique_schemes = []
    final_signatures = set()
    
    for scheme in all_schemes:
        # 生成唯一签名：核心序列 + 反向引物 + 错配位置
        signature = (
            scheme['fwd_allele1_core'],
            scheme['fwd_allele2_core'],
            scheme['reverse'],
            scheme['mismatch_pos']
        )
        
        # 如果已经见过这个签名，跳过
        if signature in final_signatures:
            continue
        
        final_signatures.add(signature)
        unique_schemes.append(scheme)
        
        # 只收集需要的数量，不强制填充
        if len(unique_schemes) >= num_schemes:
            break
    
    # 返回实际设计出的独特方案（可能少于 num_schemes）
    return unique_schemes


def design_kasp_primers_with_diagnosis(upstream: str, downstream: str, allele1: str, allele2: str,
                                        config: KASPConfig = None, num_schemes: int = 5) -> Tuple[List[Dict], Optional[DesignFailureDiagnosis]]:
    """
    KASP引物设计主入口函数 - 带智能失败诊断
    
    这是推荐使用的设计函数，会在设计失败时自动进行序列诊断并返回详细建议。
    
    ==================== 核心功能说明 ====================
    
    1. 【LGC标准人工错配 (n-3 Deliberate Mismatch)】
       - 位置：在ASP引物3'端倒数第3位引入人工错配
       - 实现：通过 apply_deliberate_mismatch() 函数
       - 规则：根据SNP强度(G/C为强, A/T为弱)选择最佳destabilizing错配
       - 作用：增强等位基因特异性，减少非特异性扩增
    
    2. 【救援模式 (Rescue Mode)】
       - 触发：标准参数无法设计出引物时自动启用
       - 参数调整：
         * 最大引物长度：25bp → 32bp
         * GC下限：35% → 20%  
         * Tm下限：55°C → 52°C
       - 实现：通过 config.get_rescue_config() 获取放宽的参数
    
    3. 【Tm贪婪延伸优化 (Tm Greedy Extension)】 - NEW!
       - 目标：复现LGC商业设计中"为了获得更好Tm而多取碱基"的行为
       - 原理：如果ASP引物Tm低于目标范围(60-62°C)，自动向5'端延伸
       - 实现：optimize_5prime_extension() 和 apply_n3_mismatch() 函数
       - 配置：TM_GREEDY_EXTENSION=True (默认启用)
       - 示例：21bp引物Tm=58°C → 延伸到22bp → Tm=60.5°C ✓
    
    4. 【智能失败诊断 (Smart Failure Feedback)】
       - 触发：救援模式也失败时
       - 分析内容：
         * GC含量分析（<25%为极低）
         * 序列复杂度评估
         * 发夹结构风险
         * 重复序列检测
       - 输出：具体失败原因和优化建议
    
    =====================================================
    
    参数:
        upstream: SNP上游序列
        downstream: SNP下游序列
        allele1, allele2: 两个等位基因碱基
        config: KASP配置参数
        num_schemes: 需要返回的方案数
    
    返回:
        (schemes, diagnosis): 
            - schemes: 引物方案列表，失败时为空列表
            - diagnosis: 失败诊断结果，成功时为None
    """
    if config is None:
        config = KASPConfig()
    
    # 调用核心设计函数
    schemes = design_kasp_primers_multi(
        upstream, downstream, allele1, allele2,
        config=config,
        num_schemes=num_schemes
    )
    
    # 如果设计成功，返回结果
    if schemes:
        return schemes, None
    
    # === 设计失败，进行智能诊断 ===
    diagnosis = diagnose_design_failure(upstream, downstream, allele1, allele2, config)
    
    return [], diagnosis


# ============================================================
# 实验指导报告生成模块 (Experimental Guide Generator)
# ============================================================

def generate_kasp_protocol(primer_data: Dict, seq_id: str = "SNP_Marker") -> str:
    """
    生成完整的KASP实验指导报告 (Markdown格式)
    
    该函数接收设计好的引物数据，返回一份可直接用于实验的详细指南，
    包含引物合成单、配方表、PCR程序等完整信息。
    
    参数:
        primer_data: 引物设计结果字典，包含以下键:
            - fwd_allele1_full: FAM标记的完整正向引物1
            - fwd_allele2_full: HEX标记的完整正向引物2
            - fwd_allele1_core: 正向引物1核心序列（不含标签）
            - fwd_allele2_core: 正向引物2核心序列（不含标签）
            - reverse: 反向引物（Common引物）
            - allele1, allele2: SNP等位基因
            - eval_fwd1, eval_fwd2, eval_rev: 引物评估信息
        seq_id: 序列/标记ID
    
    返回:
        格式化的Markdown文本，可直接用于Web渲染
    """
    
    # 提取引物信息
    fam_full = primer_data.get('fwd_allele1_full', '')
    hex_full = primer_data.get('fwd_allele2_full', '')
    fam_core = primer_data.get('fwd_allele1_core', '')
    hex_core = primer_data.get('fwd_allele2_core', '')
    common = primer_data.get('reverse', '')
    allele1 = primer_data.get('allele1', 'X')
    allele2 = primer_data.get('allele2', 'Y')
    
    # 获取Tm信息
    tm_fam = primer_data.get('eval_fwd1', {}).get('tm', 'N/A')
    tm_hex = primer_data.get('eval_fwd2', {}).get('tm', 'N/A')
    tm_common = primer_data.get('eval_rev', {}).get('tm', 'N/A')
    
    # 标准FAM/HEX尾巴（用小写标记）
    FAM_TAIL = "gaaggtgaccaagttcatgct"  # 小写表示通用标签
    HEX_TAIL = "gaaggtcggagtcaacggatt"
    
    # 构建报告
    report = f"""
# 🧪 KASP实验指导报告

**标记ID:** `{seq_id}` | **SNP类型:** `{allele1}/{allele2}` | **生成时间:** {datetime.now().strftime('%Y-%m-%d %H:%M')}

---

## 📋 1. 引物合成单 (Ready-to-Order Sequences)

> **⚠️ 请直接复制以下序列用于合成**  
> 标签部分(小写) + 引物部分(大写)，方便核对

### 🔵 Allele {allele1} - FAM标记引物 (ASP1)

```
{FAM_TAIL}{fam_core.upper()}
```

| 项目 | 信息 |
|------|------|
| 完整长度 | {len(fam_full)} bp |
| 标签长度 | 21 bp (FAM通用尾巴) |
| 核心引物 | {len(fam_core)} bp |
| Tm值 | {tm_fam}°C |

### 🟢 Allele {allele2} - HEX标记引物 (ASP2)

```
{HEX_TAIL}{hex_core.upper()}
```

| 项目 | 信息 |
|------|------|
| 完整长度 | {len(hex_full)} bp |
| 标签长度 | 21 bp (HEX通用尾巴) |
| 核心引物 | {len(hex_core)} bp |
| Tm值 | {tm_hex}°C |

### ⚪ Common引物 (通用反向引物)

```
{common.upper()}
```

| 项目 | 信息 |
|------|------|
| 长度 | {len(common)} bp |
| Tm值 | {tm_common}°C |
| 备注 | 无需加标签，直接合成 |

---

## 🧪 2. 引物混合液配方 (Primer Mix Recipe)

> 以下配方用于配制 **100 µL Assay Mix 储存液** (可进行约200个反应)

| 组分 | 体积 | 浓度 | 说明 |
|------|------|------|------|
| ASP1 (FAM-tailed) | **12 µL** | 100 µM | Allele {allele1} 特异性引物 |
| ASP2 (HEX-tailed) | **12 µL** | 100 µM | Allele {allele2} 特异性引物 |
| Common Primer | **30 µL** | 100 µM | 通用反向引物 |
| ddH₂O (无核酸酶水) | **46 µL** | - | 补足体积 |
| **总计** | **100 µL** | - | **Assay Mix** |

### 📌 配制步骤

1. 将合成的引物用 ddH₂O 稀释至 **100 µM** 储存液
2. 按上表比例混合三种引物
3. 轻轻涡旋混匀，短暂离心
4. 分装至 PCR 管中，**-20°C 长期保存**
5. 工作液可在 4°C 保存 1-2 周

---

## 🔬 3. PCR反应体系 (Reaction Setup)

> 使用 KASP Master Mix (2×)，含有 Taq酶、dNTPs、FRET Cassettes、ROX

### 单个反应体系 (10 µL 体系)

| 组分 | 体积 | 终浓度 |
|------|------|--------|
| KASP Master Mix (2×) | 5.0 µL | 1× |
| Assay Mix (上述配制) | 0.14 µL | ~12 nM ASP / ~30 nM Common |
| DNA模板 (10-20 ng/µL) | 2.0 µL | 20-40 ng |
| ddH₂O | 2.86 µL | - |
| **总计** | **10 µL** | - |

### 📌 操作提示

- **DNA浓度**: 建议使用 10-50 ng/µL 的 gDNA
- **模板量**: 每反应 10-100 ng DNA 均可
- **Master Mix**: 推荐 LGC KASP Master Mix 或兼容产品
- **高通量**: 384孔板建议使用 5 µL 体系（各组分减半）

---

## 📡 4. 仪器通道设置 (Instrument Setup)

> 适用于大多数 qPCR 仪器 (ABI, Bio-Rad, Roche等)

| 通道 | 荧光基团 | 检测目标 | 激发/发射波长 |
|------|----------|----------|---------------|
| **Channel 1** | FAM | Allele {allele1} | 495/520 nm |
| **Channel 2** | HEX/VIC | Allele {allele2} | 535/556 nm |
| **Reference** | ROX | 内参校正 (可选) | 575/602 nm |

### 📌 终点读板设置

- **读板模式**: End-point (终点读板)
- **读板温度**: 25-37°C
- **读板时间**: 扩增完成后静置 1 min 再读板

---

## 🌡️ 5. 降落PCR程序 (Touchdown Protocol)

> **标准KASP循环程序** - 适用于大多数SNP标记

```
┌─────────────────────────────────────────────────────────┐
│  阶段1: 预变性 (热启动酶激活)                              │
│  ─────────────────────────────────                      │
│  94°C ─── 15 min ─── (激活热启动Taq酶)                   │
├─────────────────────────────────────────────────────────┤
│  阶段2: 降落扩增 (Touchdown) × 10 循环                    │
│  ─────────────────────────────────                      │
│  94°C ─── 20 sec ─── (变性)                              │
│  61°C ─── 60 sec ─── (退火/延伸，每循环 -0.6°C)           │
│                                                         │
│  📉 退火温度变化: 61→60.4→59.8→...→55.6°C                │
├─────────────────────────────────────────────────────────┤
│  阶段3: 常规扩增 × 26-35 循环                             │
│  ─────────────────────────────────                      │
│  94°C ─── 20 sec ─── (变性)                              │
│  55°C ─── 60 sec ─── (退火/延伸) ⭐ 荧光信号采集          │
├─────────────────────────────────────────────────────────┤
│  阶段4: 终点读板                                          │
│  ─────────────────────────────────                      │
│  37°C ─── 1 min ─── (End-point Read)                    │
│                                                         │
│  📊 此时采集最终荧光数据用于基因分型                        │
└─────────────────────────────────────────────────────────┘
```

### 📌 程序参数速查表

| 阶段 | 温度 | 时间 | 循环数 | 说明 |
|------|------|------|--------|------|
| 预变性 | 94°C | 15 min | 1× | 激活热启动酶 |
| 降落扩增 | 94°C | 20 sec | 10× | 变性 |
| | 61°C→55°C | 60 sec | (-0.6°C/循环) | 降落退火 |
| 常规扩增 | 94°C | 20 sec | 26-35× | 变性 |
| | 55°C | 60 sec | ⭐采集荧光 | 退火/延伸 |
| 终点读板 | 37°C | 1 min | 1× | End-point |

---

## 📊 6. 结果判读 (Result Interpretation)

### 典型散点图分布

```
       HEX荧光 (Allele {allele2})
          ↑
          │    ● ● ●          (纯合 {allele2}/{allele2})
          │      ● ●
          │
          │         ● ● ●     (杂合 {allele1}/{allele2})
          │           ● ●
          │
          │                   (纯合 {allele1}/{allele1})
          │              ● ● ●
          │                ● ●
          └──────────────────────→ FAM荧光 (Allele {allele1})
```

### 分型标准

| 位置 | 基因型 | 说明 |
|------|--------|------|
| **右下角** (高FAM/低HEX) | `{allele1}/{allele1}` | 纯合型1 |
| **左上角** (低FAM/高HEX) | `{allele2}/{allele2}` | 纯合型2 |
| **中间** (中FAM/中HEX) | `{allele1}/{allele2}` | 杂合型 |
| **左下角** (低FAM/低HEX) | NTC/失败 | 无模板对照或扩增失败 |

---

## ⚠️ 7. 常见问题排查

| 问题 | 可能原因 | 解决方案 |
|------|----------|----------|
| 无荧光信号 | DNA量不足/降解 | 检查DNA浓度和完整性 |
| 所有样品聚成一团 | 退火温度过高/低 | 调整Touchdown起始温度 |
| 杂合型分离不清 | 循环数不足 | 增加5-10个循环 |
| 背景荧光高 | 引物浓度过高 | 减少Assay Mix用量 |
| 批次间差异大 | 引物或Master Mix批次不同 | 同一实验使用同批次试剂 |

---

## 📦 8. 试剂采购清单

| 试剂 | 推荐品牌/型号 | 用量估算 |
|------|---------------|----------|
| KASP Master Mix (2×) | LGC Biosearch / KBS-1016-002 | 1 mL = 200反应 |
| 引物合成 | 华大基因 / 生工生物 / Invitrogen | 按上述序列 |
| 96/384孔板 | Applied Biosystems / Bio-Rad | PCR专用白色板 |
| 封板膜 | 透明光学封板膜 | 适配qPCR仪 |

---

*📋 本报告由 KASP引物设计工具 v7.0 自动生成 | 基于LGC标准KASP流程*
"""
    
    return report


def generate_kasp_protocol_brief(primer_data: Dict, seq_id: str = "SNP_Marker") -> str:
    """
    生成简化版KASP实验指导 (用于快速参考)
    
    参数:
        primer_data: 引物设计结果字典
        seq_id: 序列ID
    
    返回:
        简化版Markdown文本
    """
    fam_full = primer_data.get('fwd_allele1_full', '')
    hex_full = primer_data.get('fwd_allele2_full', '')
    common = primer_data.get('reverse', '')
    allele1 = primer_data.get('allele1', 'X')
    allele2 = primer_data.get('allele2', 'Y')
    
    # 分离标签和核心序列（用大小写区分）
    FAM_TAIL = "gaaggtgaccaagttcatgct"
    HEX_TAIL = "gaaggtcggagtcaacggatt"
    fam_core = primer_data.get('fwd_allele1_core', '').upper()
    hex_core = primer_data.get('fwd_allele2_core', '').upper()
    
    brief = f"""
### 🧬 引物合成单 - {seq_id}

**直接复制以下序列用于合成** (小写=通用标签, 大写=特异性引物)

| 引物名称 | 序列 (5'→3') | 长度 |
|----------|--------------|------|
| **{seq_id}_FAM** ({allele1}) | `{FAM_TAIL}{fam_core}` | {len(fam_full)} bp |
| **{seq_id}_HEX** ({allele2}) | `{HEX_TAIL}{hex_core}` | {len(hex_full)} bp |
| **{seq_id}_Common** | `{common.upper()}` | {len(common)} bp |

---

### 🧪 快速配方 (100 µL Assay Mix)

```
ASP1 (FAM):    12 µL  (100 µM)
ASP2 (HEX):    12 µL  (100 µM)  
Common:        30 µL  (100 µM)
ddH₂O:         46 µL
────────────────────────────
Total:        100 µL
```

### ⏱️ PCR程序速记

```
94°C/15min → [94°C/20s + 61→55°C/60s]×10 → [94°C/20s + 55°C/60s]×30 → 37°C/1min
```
"""
    return brief


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
    """严格评估引物质量（用于常规PCR）- 使用Primer3热力学分析"""
    seq = seq.upper()
    seq = re.sub(r'[^ATGC]', '', seq)
    
    gc_content = calc_gc_content(seq)
    tm = calc_tm(seq)
    
    # 使用Primer3进行结构分析
    hairpin_result = check_hairpin(seq)
    homodimer_result = check_homodimer(seq)
    
    has_hairpin = hairpin_result[0] if isinstance(hairpin_result, tuple) else hairpin_result
    hairpin_tm = hairpin_result[1] if isinstance(hairpin_result, tuple) and len(hairpin_result) > 1 else None
    
    has_self_dimer = homodimer_result[0] if isinstance(homodimer_result, tuple) else homodimer_result
    homodimer_dg = homodimer_result[1] if isinstance(homodimer_result, tuple) and len(homodimer_result) > 1 else None
    
    three_prime_ok, three_prime_msg = check_gc_clamp(seq)
    has_repeat = check_repeat_region(seq)
    
    result = {
        'sequence': seq,
        'length': len(seq),
        'gc_content': gc_content,
        'tm': tm,
        'has_hairpin': has_hairpin,
        'hairpin_tm': hairpin_tm,
        'has_self_dimer': has_self_dimer,
        'homodimer_dg': homodimer_dg,
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


# ==================== 基于Primer3的引物设计 (polyoligo-kasp风格) ====================

def design_primers_with_primer3(sequence: str, config: RegularPCRConfig = None,
                                  num_pairs: int = 5, target_region: Tuple[int, int] = None) -> List[Dict]:
    """
    使用Primer3库进行专业的引物设计 - 参考polyoligo-kasp设计理念
    
    参数:
        sequence: 模板序列
        config: PCR配置参数
        num_pairs: 需要返回的引物对数量
        target_region: 目标区域 (start, length)
    
    返回:
        引物对列表
    """
    if not PRIMER3_AVAILABLE:
        return []  # 回退到手动设计
    
    if config is None:
        config = RegularPCRConfig()
    
    sequence = re.sub(r'[^ATGC]', '', sequence.upper())
    seq_len = len(sequence)
    
    if seq_len < config.PRODUCT_MIN + 40:
        return []
    
    # 构建Primer3序列参数
    seq_args = {
        'SEQUENCE_ID': 'target',
        'SEQUENCE_TEMPLATE': sequence,
    }
    
    # 设置目标区域
    if target_region:
        seq_args['SEQUENCE_TARGET'] = [target_region[0], target_region[1]]
    else:
        # 默认目标区域为序列中间
        margin = max(50, config.MIN_PRIMER_LEN + 10)
        target_start = margin
        target_len = seq_len - 2 * margin
        if target_len > 0:
            seq_args['SEQUENCE_INCLUDED_REGION'] = [target_start, target_len]
    
    # 构建Primer3全局参数 - 参考polyoligo-kasp
    global_args = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_NUM_RETURN': num_pairs * 2,  # 多返回一些以便筛选
        
        # 引物长度
        'PRIMER_MIN_SIZE': config.MIN_PRIMER_LEN,
        'PRIMER_OPT_SIZE': config.OPTIMAL_PRIMER_LEN,
        'PRIMER_MAX_SIZE': config.MAX_PRIMER_LEN,
        
        # Tm参数
        'PRIMER_MIN_TM': config.MIN_TM,
        'PRIMER_OPT_TM': config.OPTIMAL_TM,
        'PRIMER_MAX_TM': config.MAX_TM,
        'PRIMER_PAIR_MAX_DIFF_TM': config.MAX_TM_DIFF,
        
        # GC含量
        'PRIMER_MIN_GC': config.MIN_GC,
        'PRIMER_MAX_GC': config.MAX_GC,
        
        # 产物大小范围
        'PRIMER_PRODUCT_SIZE_RANGE': [[config.PRODUCT_MIN, config.PRODUCT_MAX]],
        
        # 热力学参数
        'PRIMER_SALT_MONOVALENT': THERMO_PARAMS.mv_conc,
        'PRIMER_SALT_DIVALENT': THERMO_PARAMS.dv_conc,
        'PRIMER_DNTP_CONC': THERMO_PARAMS.dntp_conc,
        'PRIMER_DNA_CONC': THERMO_PARAMS.dna_conc,
        
        # 二级结构和二聚体检测阈值
        'PRIMER_MAX_SELF_ANY': 8,
        'PRIMER_MAX_SELF_END': 3,
        'PRIMER_PAIR_MAX_COMPL_ANY': 8,
        'PRIMER_PAIR_MAX_COMPL_END': 3,
        'PRIMER_MAX_HAIRPIN_TH': 47.0,
        
        # 其他约束
        'PRIMER_MAX_POLY_X': 4,  # 最多4个连续相同碱基
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_GC_CLAMP': 1,  # 3'端至少1个GC
        
        # Tm计算方法
        'PRIMER_TM_FORMULA': 1,  # SantaLucia
        'PRIMER_SALT_CORRECTIONS': 1,  # SantaLucia
    }
    
    try:
        # 调用Primer3设计
        results = primer3.design_primers(seq_args, global_args)
        
        # 解析结果
        primer_pairs = []
        num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
        
        for i in range(min(num_returned, num_pairs)):
            try:
                # 提取引物信息
                left_seq = results.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
                right_seq = results.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
                
                if not left_seq or not right_seq:
                    continue
                
                left_pos = results.get(f'PRIMER_LEFT_{i}', [0, 0])
                right_pos = results.get(f'PRIMER_RIGHT_{i}', [0, 0])
                
                left_tm = results.get(f'PRIMER_LEFT_{i}_TM', 0)
                right_tm = results.get(f'PRIMER_RIGHT_{i}_TM', 0)
                
                left_gc = results.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0)
                right_gc = results.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0)
                
                product_size = results.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0)
                
                # 检测二聚体
                has_dimer, dimer_dg, dimer_tm = check_heterodimer(left_seq, right_seq)
                
                # 检测发夹
                left_hairpin, left_hp_tm = check_hairpin(left_seq)
                right_hairpin, right_hp_tm = check_hairpin(right_seq)
                
                # 检测自身二聚体
                left_homo, left_homo_dg, _ = check_homodimer(left_seq)
                right_homo, right_homo_dg, _ = check_homodimer(right_seq)
                
                # 构建评估结果
                fwd_eval = {
                    'sequence': left_seq,
                    'length': len(left_seq),
                    'tm': round(left_tm, 1),
                    'gc_content': round(left_gc, 1),
                    'has_hairpin': left_hairpin,
                    'has_self_dimer': left_homo,
                    'hairpin_tm': left_hp_tm,
                    'homodimer_dg': left_homo_dg,
                    'issues': [],
                    'score': 100
                }
                
                rev_eval = {
                    'sequence': right_seq,
                    'length': len(right_seq),
                    'tm': round(right_tm, 1),
                    'gc_content': round(right_gc, 1),
                    'has_hairpin': right_hairpin,
                    'has_self_dimer': right_homo,
                    'hairpin_tm': right_hp_tm,
                    'homodimer_dg': right_homo_dg,
                    'issues': [],
                    'score': 100
                }
                
                # 评分调整
                for eval_result in [fwd_eval, rev_eval]:
                    if eval_result['has_hairpin']:
                        eval_result['issues'].append('可能形成发夹结构')
                        eval_result['score'] -= 10
                    if eval_result['has_self_dimer']:
                        eval_result['issues'].append('自身二聚体风险')
                        eval_result['score'] -= 10
                
                # 计算综合评分
                tm_diff = abs(left_tm - right_tm)
                base_score = (fwd_eval['score'] + rev_eval['score']) / 2
                
                if tm_diff <= 0.5:
                    base_score += 10
                elif tm_diff <= 1.0:
                    base_score += 5
                elif tm_diff > 2.0:
                    base_score -= 10
                
                if has_dimer:
                    base_score -= 15
                
                # 产物大小评分
                if 200 <= product_size <= 400:
                    base_score += 5
                
                total_score = max(0, min(100, base_score))
                
                # 可用性判断
                is_usable = (
                    not has_dimer and
                    tm_diff <= config.MAX_TM_DIFF and
                    config.MIN_TM <= left_tm <= config.MAX_TM and
                    config.MIN_TM <= right_tm <= config.MAX_TM
                )
                
                pair = {
                    'forward': left_seq,
                    'reverse': right_seq,
                    'fwd_start': left_pos[0] + 1,
                    'fwd_end': left_pos[0] + left_pos[1],
                    'rev_start': right_pos[0] - right_pos[1] + 2,
                    'rev_end': right_pos[0] + 1,
                    'fwd_eval': fwd_eval,
                    'rev_eval': rev_eval,
                    'tm_diff': round(tm_diff, 1),
                    'has_dimer': has_dimer,
                    'heterodimer_dg': dimer_dg,
                    'product_size': product_size,
                    'total_score': total_score,
                    'is_usable': is_usable,
                    'wheat_mode': config.WHEAT_MODE,
                    'wheat_issues': [],
                    'primer3_penalty': results.get(f'PRIMER_PAIR_{i}_PENALTY', 0),
                    'fwd_position_percent': round(left_pos[0] / seq_len * 100, 1),
                    'rev_position_percent': round(right_pos[0] / seq_len * 100, 1)
                }
                
                primer_pairs.append(pair)
                
            except Exception:
                continue
        
        # 按评分排序
        primer_pairs.sort(key=lambda x: x['total_score'], reverse=True)
        
        return primer_pairs[:num_pairs]
        
    except Exception as e:
        return []


def design_regular_primers(sequence: str, config: RegularPCRConfig = None, 
                          num_pairs: int = 5, target_start: int = None, 
                          target_end: int = None) -> List[Dict]:
    """
    设计常规PCR引物对 - 优化版（支持小麦模式）
    优先使用Primer3库，失败时回退到手动设计
    确保不产生重复引物对，质量不达标时返回空列表
    """
    if config is None:
        config = RegularPCRConfig()
    
    sequence = re.sub(r'[^ATGC]', '', sequence.upper())
    seq_len = len(sequence)
    
    # 序列太短直接返回空
    min_required_len = config.PRODUCT_MIN + 2 * config.MIN_PRIMER_LEN
    if seq_len < min_required_len:
        return []  # 序列太短，无法设计
    
    # === 优先尝试Primer3设计 (非小麦模式或小麦模式均可) ===
    if PRIMER3_AVAILABLE:
        target_region = None
        if target_start is not None and target_end is not None:
            target_len = target_end - target_start
            if target_len > 10:
                target_region = (target_start, target_len)
        
        primer3_results = design_primers_with_primer3(
            sequence, config, num_pairs, target_region
        )
        
        # 如果Primer3返回了足够的结果，进行小麦特异性评估后返回
        if primer3_results and len(primer3_results) >= max(1, num_pairs // 2):
            # 对小麦模式进行额外评估
            if config.WHEAT_MODE:
                for pair in primer3_results:
                    fwd_position = pair['fwd_start'] / seq_len
                    rev_position = pair['rev_end'] / seq_len
                    
                    fwd_bonus, fwd_issues = check_wheat_specificity(pair['forward'], fwd_position)
                    rev_bonus, rev_issues = check_wheat_specificity(pair['reverse'], rev_position)
                    
                    pair['wheat_issues'] = fwd_issues + rev_issues
                    pair['total_score'] = max(0, min(100, pair['total_score'] + (fwd_bonus + rev_bonus) / 2))
                    
                    # 两个引物都在3'半区给额外奖励
                    if fwd_position > 0.5 and rev_position > 0.6:
                        pair['total_score'] = min(100, pair['total_score'] + 10)
                
                # 重新排序
                primer3_results.sort(key=lambda x: x['total_score'], reverse=True)
            
            return primer3_results
    
    # === 回退到手动设计 ===
    
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
    
    # 确保目标区域有效
    if target_end <= target_start or target_end - target_start < config.PRODUCT_MIN:
        return []  # 目标区域无效
    
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
    """生成KASP引物CSV内容 - 增强版（包含人工错配和救援模式信息）"""
    output = io.StringIO()
    writer = csv.writer(output)
    
    writer.writerow(['KASP引物设计报告 (优化版 v7.0)'])
    writer.writerow(['序列ID', seq_id])
    writer.writerow(['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
    
    # 检查是否为小麦模式
    is_wheat_mode = schemes[0].get('wheat_mode', False) if schemes else False
    has_rescue = any(s.get('rescue_mode', False) for s in schemes)
    
    if is_wheat_mode:
        writer.writerow(['模式', '🌾 小麦KASP模式 (五大忌检测)'])
    if has_rescue:
        writer.writerow(['提示', '⚠️ 部分引物由救援模式设计（用于AT-rich序列）'])
    writer.writerow(['人工错配', 'LGC标准 - 在n-3位置引入Deliberate Mismatch增强特异性'])
    writer.writerow([])
    
    # 表头
    headers = ['方案', '评分', '等级', '可用性', '设计模式',
               'FAM引物(完整)', 'HEX引物(完整)', '通用反向引物',
               'Allele1-Tm', 'Allele2-Tm', 'ASP Tm差', 'ASP-Common Tm差',
               'Allele1-GC%', 'Allele2-GC%', 'Rev-Tm', 'Rev-GC%',
               '产物大小', '错配位置', '错配变化', 'SNP类型']
    
    if is_wheat_mode:
        headers.extend(['小麦评估', '注意事项'])
    
    writer.writerow(headers)
    
    for i, scheme in enumerate(schemes, 1):
        grade, stars, _ = get_quality_grade(scheme['total_score'])
        is_usable = scheme.get('is_usable', True)
        rescue_mode = scheme.get('rescue_mode', False)
        
        # 获取人工错配详细信息
        mismatch_info = scheme.get('deliberate_mismatch_info', {})
        snp_bases = mismatch_info.get('snp_bases', (scheme.get('allele1', '?'), scheme.get('allele2', '?')))
        snp_type = f"[{snp_bases[0]}/{snp_bases[1]}]"
        
        asp_common_diff = scheme.get('asp_common_tm_diff', 0)
        
        row = [
            f"方案{i}", f"{scheme['total_score']:.1f}", f"{grade} {stars}",
            "推荐" if is_usable else "慎用",
            "救援模式" if rescue_mode else "标准模式",
            scheme['fwd_allele1_full'], scheme['fwd_allele2_full'], scheme['reverse'],
            f"{scheme['eval_fwd1']['tm']}°C", f"{scheme['eval_fwd2']['tm']}°C", 
            f"{scheme['tm_diff']:.1f}°C", f"{asp_common_diff:.1f}°C" if asp_common_diff else "N/A",
            f"{scheme['eval_fwd1']['gc_content']:.1f}%", f"{scheme['eval_fwd2']['gc_content']:.1f}%",
            f"{scheme['eval_rev']['tm']}°C", f"{scheme['eval_rev']['gc_content']:.1f}%",
            f"{scheme['product_size']}bp", f"n{scheme['mismatch_pos']}", scheme['mismatch_change'],
            snp_type
        ]
        
        if is_wheat_mode:
            amplicon_status = scheme.get('wheat_details', {}).get('amplicon_status', '')
            wheat_issues = scheme.get('wheat_issues', [])
            row.append(amplicon_status)
            row.append("; ".join(wheat_issues) if wheat_issues else "无")
        
        writer.writerow(row)
    
    # 人工错配说明
    writer.writerow([])
    writer.writerow(['=== LGC人工错配(Deliberate Mismatch)说明 ==='])
    writer.writerow(['位置', 'n-3位置', '在ASP引物3\'端倒数第3位引入错配'])
    writer.writerow(['作用', '增强特异性', '使非目标等位基因的引物结合更不稳定'])
    writer.writerow(['规则', '强SNP(G/C)', '使用强destabilizing错配(如A→C)'])
    writer.writerow(['规则', '弱SNP(A/T)', '使用中等destabilizing错配(如A→G)'])
    
    # 救援模式说明
    if has_rescue:
        writer.writerow([])
        writer.writerow(['=== 救援模式说明 ==='])
        writer.writerow(['触发条件', 'AT-rich序列', '标准参数无法设计出合格引物时自动触发'])
        writer.writerow(['参数调整', '引物长度', '最大延长至30-32bp'])
        writer.writerow(['参数调整', 'GC下限', '降低至20%'])
        writer.writerow(['参数调整', 'Tm下限', '降低至52°C'])
        writer.writerow(['注意', '实验验证', '救援模式引物需更多实验验证'])
    
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
    
    # 初始化会话状态
    if 'last_snp_input' not in st.session_state:
        st.session_state['last_snp_input'] = None
    if 'restore_last' not in st.session_state:
        st.session_state['restore_last'] = False
    
    # 初始化临时编辑缓冲区（用于页面切换时保存未提交的输入）
    if 'temp_seq_input' not in st.session_state:
        st.session_state['temp_seq_input'] = None
    if 'temp_seq_id' not in st.session_state:
        st.session_state['temp_seq_id'] = None
    if 'temp_num_schemes' not in st.session_state:
        st.session_state['temp_num_schemes'] = 5
    if 'temp_wheat_mode' not in st.session_state:
        st.session_state['temp_wheat_mode'] = False
    if 'temp_min_primer_len' not in st.session_state:
        st.session_state['temp_min_primer_len'] = 18
    if 'temp_max_primer_len' not in st.session_state:
        st.session_state['temp_max_primer_len'] = 30
    if 'temp_min_tm' not in st.session_state:
        st.session_state['temp_min_tm'] = 55.0
    if 'temp_max_tm' not in st.session_state:
        st.session_state['temp_max_tm'] = 68.0
    
    # 检查是否有保存的上次输入
    has_saved_input = st.session_state['last_snp_input'] is not None
    
    # 显示临时输入恢复提示（如果有临时保存的输入且不是从上次输入恢复）
    has_temp_input = st.session_state['temp_seq_input'] is not None and not st.session_state.get('restore_last', False)
    if has_temp_input:
        st.info(f"💫 检测到未保存的编辑内容 | [{st.session_state['temp_seq_id'] or '未命名'}]")
        temp_col1, temp_col2 = st.columns([1, 1])
        with temp_col1:
            if st.button("✏️ 恢复编辑", use_container_width=True, key="restore_temp_input"):
                st.session_state['restore_last'] = False
        with temp_col2:
            if st.button("🗑️ 清除草稿", use_container_width=True, key="clear_temp_input"):
                st.session_state['temp_seq_input'] = None
                st.session_state['temp_seq_id'] = None
                st.rerun()
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # 显示恢复选项
        if has_saved_input:
            saved_time = st.session_state['last_snp_input'].get('timestamp', '未知时间')
            saved_id = st.session_state['last_snp_input'].get('seq_id', '未命名')
            restore_col1, restore_col2, restore_col3 = st.columns([2.5, 0.8, 0.7])
            with restore_col1:
                st.info(f"💾 上次输入: {saved_id} ({saved_time})")
            with restore_col2:
                if st.button("📂 恢复", help="恢复上次输入的序列和参数", key="restore_last_input", use_container_width=True):
                    st.session_state['restore_last'] = True
            with restore_col3:
                if st.button("🗑️ 清除", help="清除保存的输入", key="clear_last_input", use_container_width=True):
                    st.session_state['last_snp_input'] = None
                    st.rerun()
        
        # 判断使用恢复的值、临时保存值还是示例
        if has_saved_input and st.session_state.get('restore_last', False):
            default_seq = st.session_state['last_snp_input'].get('sequence', example_seq)
            default_id = st.session_state['last_snp_input'].get('seq_id', 'My_SNP_Marker')
        elif st.session_state['temp_seq_input'] is not None:
            default_seq = st.session_state['temp_seq_input']
            default_id = st.session_state['temp_seq_id'] or 'My_SNP_Marker'
        else:
            default_seq = example_seq
            default_id = 'My_SNP_Marker'
        
        seq_input = st.text_area(
            "输入序列（包含SNP标记）",
            value=default_seq,
            height=200,
            help="SNP位点使用 [碱基1/碱基2] 格式标记",
            key="seq_input_field"
        )
        # 自动保存当前编辑内容到临时缓冲区
        st.session_state['temp_seq_input'] = seq_input
        
        seq_id = st.text_input("序列ID（可选）", value=default_id, key="seq_id_field")
        # 自动保存序列ID到临时缓冲区
        st.session_state['temp_seq_id'] = seq_id
    
    with col2:
        st.markdown("**参数设置**")
        num_schemes = st.slider("生成方案数", 3, 15, st.session_state.get('temp_num_schemes', 5), key="num_schemes_slider")
        # 保存到临时缓冲区
        st.session_state['temp_num_schemes'] = num_schemes
        
        # 🌾 小麦模式
        st.markdown("---")
        wheat_mode = st.checkbox("🌾 小麦KASP模式", value=st.session_state.get('temp_wheat_mode', False),
                                  help="针对小麦六倍体(AABBDD)优化，检测五大忌", key="wheat_mode_check")
        # 保存到临时缓冲区
        st.session_state['temp_wheat_mode'] = wheat_mode
        
        if wheat_mode:
            st.warning("""**🌾 小麦五大忌检测已启用：**
1️⃣ 同源基因干扰 → 请BLAST验证
2️⃣ 侧翼SNP干扰 → 自动检测
3️⃣ 扩增子过长 → 限制50-100bp
4️⃣ GC含量极端 → 检测30-65%
5️⃣ 重复序列 → 检测转座子/SSR""")
        
        with st.expander("高级参数"):
            min_primer_len = st.number_input("最小引物长度", 15, 25, st.session_state.get('temp_min_primer_len', 18), key="min_primer_len_input")
            st.session_state['temp_min_primer_len'] = min_primer_len
            max_primer_len = st.number_input("最大引物长度", 20, 35, st.session_state.get('temp_max_primer_len', 30), key="max_primer_len_input")
            st.session_state['temp_max_primer_len'] = max_primer_len
            min_tm = st.number_input("最低Tm (°C)", 50.0, 60.0, st.session_state.get('temp_min_tm', 55.0), key="min_tm_input")
            st.session_state['temp_min_tm'] = min_tm
            max_tm = st.number_input("最高Tm (°C)", 60.0, 75.0, st.session_state.get('temp_max_tm', 68.0), key="max_tm_input")
            st.session_state['temp_max_tm'] = max_tm
            
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
                
                # === 使用带智能诊断的设计函数 ===
                # 核心设计流程说明:
                # 1. 首先尝试标准模式设计（含n-3人工错配）
                # 2. 如果失败，自动启用救援模式（放宽参数）
                # 3. 如果救援模式也失败，进行智能诊断
                schemes, diagnosis = design_kasp_primers_with_diagnosis(
                    upstream, downstream, allele1, allele2, config, num_schemes
                )
            
            # === 智能失败诊断显示 ===
            if not schemes:
                st.error("❌ 未能设计出合适的引物")
                
                # 显示详细诊断信息
                if diagnosis:
                    # 根据严重程度选择颜色
                    if diagnosis.severity == 'critical':
                        st.error(f"🔴 **主要原因**: {diagnosis.failure_reason}")
                    elif diagnosis.severity == 'high':
                        st.warning(f"🟠 **主要原因**: {diagnosis.failure_reason}")
                    else:
                        st.info(f"🟡 **主要原因**: {diagnosis.failure_reason}")
                    
                    # 显示格式化的诊断报告
                    st.markdown(format_diagnosis_for_display(diagnosis))
                    
                    # 可折叠的详细数据
                    with st.expander("📋 查看原始诊断数据"):
                        st.json(diagnosis.sequence_analysis)
                else:
                    st.warning("未能获取详细诊断信息，请检查输入序列格式")
                
                return
            
            # 显示实际设计出的方案数量
            actual_count = len(schemes)
            if actual_count < num_schemes:
                st.success(f"✅ 成功设计 {actual_count} 套独特引物方案（已去除重复方案）")
                st.info(f"ℹ️ 由于序列特性限制，实际可设计的独特方案少于请求的 {num_schemes} 套")
            else:
                st.success(f"✅ 成功设计 {actual_count} 套引物方案！")
            
            # 小麦模式警告
            if wheat_mode:
                st.info("""**🌾 重要提醒 - 大忌#1 同源基因干扰：**
请将设计的引物序列BLAST到小麦A、B、D三个基因组，确认：
- 引物是否只在目标基因组有完美匹配
- 是否需要在Common引物区域使用Genome-specific SNP
- 推荐工具：[PolyMarker](http://polymarker.tgac.ac.uk/) | [CerealsDB](http://www.cerealsdb.uk.net/)""")
            
            # 生物学验证提示 (非小麦模式也显示)
            st.markdown("---")
            st.markdown("### 🧪 生物学验证（确保引物特异性）")
            
            # 为小麦模式特殊显示验证要求
            if wheat_mode:
                st.warning("""
## 🌾 小麦KASP引物 - 必须验证！

你正在使用**小麦KASP模式**。小麦是**六倍体 (AABBDD)**，具有高度同源的A、B、D三个基因组。
**不经过BLAST验证的引物设计本质上是不完整的**。

**为什么必须验证？** 虽然本工具已经进行了初步的同源检测，但最终确认还需要在实际的参考基因组上验证。
                """)
            
            with st.expander("📍 Ensembl Plants BLAST验证步骤（强烈推荐）", expanded=wheat_mode):
                st.markdown("""
#### 第一步：打开BLAST工具

选择一个工具进行验证（推荐顺序）：

1. **🏆 PolyMarker（最佳选择，专为小麦设计）**
   - 地址：http://polymarker.tgac.ac.uk/
   - 优点：自动检查小麦A/B/D三个基因组，一次性出结果
   - 推荐用于：小麦KASP引物验证

2. **Ensembl Plants BLAST（官方工具）**
   - 地址：https://plants.ensembl.org/Multi/Tools/Blast?db=core
   - 优点：官方权威，支持多个物种和基因组版本
   - 缺点：需要手动检查每个基因组

3. **NCBI BLAST（通用工具）**
   - 地址：https://blast.ncbi.nlm.nih.gov/
   - 优点：数据库全面
   - 缺点：需要选择正确的数据库

---

#### 第二步：输入序列

- 复制设计好的 **反向引物（Common引物）序列** 或 **PCR产物序列**
- 粘贴到BLAST的输入框中
- 如果用PolyMarker，粘贴反向引物序列即可

---

#### 第三步：设置BLAST参数（以Ensembl为例）

**对于小麦序列：**
```
Search against: Wheat (Triticum aestivum)
Database: IWGSC (v2.1)        ← 这是最新的小麦参考基因组
Program: megablast             ← 用于相同或高度相似的序列
Other settings: Keep default   ← 默认参数通常适合
```

**关键要点：**
- 确保选择的是 **IWGSC (v2.1)** 而不是其他版本
- 对于 **Allele特异性引物**，单独验证该引物序列
- 对于 **Common引物**，验证该反向引物序列

---

#### 第四步：运行BLAST

- 点击 "Run BLAST" 或 "Search" 按钮
- 等待搜索完成（通常几秒到几分钟）

---

#### 第五步：分析结果（这一步决定生死！）

**🎯 对于小麦（使用PolyMarker或分别检查A/B/D基因组）：**

✅ **理想情况（可以用）：**
```
目标染色体 2B：
  ├─ 100% Identity, 0 Gaps          ← 完美匹配 ✓
  └─ Length: 20bp (比如你的引物长度)

其他染色体：
  ├─ 2A: <90% Identity 或有Gap      ← 不存在干扰
  ├─ 2D: <90% Identity 或有Gap      ← 不存在干扰
  └─ 其他: 无匹配或很低匹配          ← 特异性良好 ✓✓✓
```

⚠️ **警告情况（需要修改引物）：**
```
2A 或 2D 上看到 99%-100% 匹配
  ↓
说明这条引物会同时扩增多个基因组
  ↓
结果：无法区分基因型，设计失败 ✗
  ↓
解决方案：
  1. 在引物3'端引入SNP（强制不匹配其他基因组）
  2. 或者选择Genome-specific SNP位点
  3. 或重新选择引物设计方案
```

❌ **不可接受的情况（必须放弃）：**
- 有多条高匹配（99%-100%），分散在不同位置
- 中间出现大的Gap，表示序列可能在重复区域
- 无法清晰识别目标匹配位置

---

#### 第六步（可选）：使用PolyMarker深度验证

如果使用PolyMarker，输出界面会直接显示：
```
引物名称: Your_Primer
Matches:
  2B_v2.1:     100% (Perfect match) ✓
  2A_v2.1:     85%  (Not matching) ✓
  2D_v2.1:     82%  (Not matching) ✓
  
Conclusion: SUITABLE FOR KASP ✅
```

这种情况下，你的引物设计合格 ✓✓✓

---

### 📋 验证清单

在使用引物前，请确认：

- [ ] 通过了 Ensembl BLAST 或 PolyMarker 验证
- [ ] 目标基因组上 100% 完美匹配
- [ ] 其他基因组上 <90% 匹配或无匹配
- [ ] 对 Allele 特异性引物也进行了验证
- [ ] 记录了 BLAST 结果（以备后续查证）

---

### ⚡ 常见问题

**Q: PolyMarker和Ensembl哪个更准确？**
A: PolyMarker专为小麦设计，用户界面更友好，推荐首选。Ensembl是官方工具，数据更新可能稍晚。

**Q: 如果两条引物验证结果不一致怎么办？**
A: 这很正常。两条引物需要分别验证。只要两条都通过验证（各自在正确位置有完美匹配），就可以组合使用。

**Q: 引物在多个位置有匹配怎么办？**
A: 如果其他位置的匹配度 <90%，仍然可以接受。但如果有多个 99%-100% 的匹配，说明引物特异性有问题。

**Q: IWGSC版本是最新的吗？**
A: v2.1 是当前最新的小麦参考基因组。如果有新版本发布，请更新。

---

**重要提醒：** 本工具的引物设计质量很大程度上取决于此验证步骤。
不经过BLAST验证的引物，即使在本工具中评分再高，在实验中也可能失败。
                """)
            
            # 保存输入信息到会话状态
            st.session_state['last_snp_input'] = {
                'sequence': seq_input,
                'seq_id': seq_id,
                'wheat_mode': wheat_mode,
                'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'upstream': upstream,
                'downstream': downstream,
                'allele1': allele1,
                'allele2': allele2
            }
            
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
                    
                    # 救援模式提示
                    if scheme.get('rescue_mode', False):
                        st.warning("🆘 **救援模式引物**: 该引物由放宽参数设计（用于AT-rich序列），建议优先验证实验效果")
                    
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
                        # ASP间Tm差异
                        tm_status = "✓" if scheme['tm_diff'] <= 1.5 else ("△" if scheme['tm_diff'] <= 2.0 else "✗")
                        st.write(f"ASP Tm差异: {scheme['tm_diff']:.1f}°C {tm_status}")
                        
                        # ASP-Common Tm差异（新增）
                        asp_common_diff = scheme.get('asp_common_tm_diff', 0)
                        if asp_common_diff:
                            tm_balance_status = "✓" if asp_common_diff <= 2.0 else ("△" if asp_common_diff <= 3.0 else "✗")
                            st.write(f"ASP-Common Tm差: {asp_common_diff:.1f}°C {tm_balance_status}")
                        
                        # 人工错配信息（增强版）
                        mismatch_info = scheme.get('deliberate_mismatch_info', {})
                        if mismatch_info:
                            snp_bases = mismatch_info.get('snp_bases', ('?', '?'))
                            st.write(f"🔬 **人工错配 (LGC标准)**")
                            st.write(f"  位置: n{scheme['mismatch_pos']} | SNP: [{snp_bases[0]}/{snp_bases[1]}]")
                            st.write(f"  变化: {scheme['mismatch_change']}")
                        else:
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
            
            # ============================================================
            # 🧪 实验指导报告 (Experimental Guide)
            # ============================================================
            st.markdown("---")
            st.markdown("## 🧪 实验指导报告")
            st.info("📘 以下是基于LGC标准KASP流程的完整实验指南，适合初学者参考使用")
            
            # 生成实验指导报告
            protocol_report = generate_kasp_protocol(best, seq_id if seq_id else "SNP_Marker")
            
            # 使用Tab展示不同版本的报告
            tab_full, tab_brief = st.tabs(["📖 完整指南", "📋 快速参考"])
            
            with tab_full:
                st.markdown(protocol_report)
                
                # 提供下载完整报告的按钮
                st.download_button(
                    label="📥 下载完整实验指南 (Markdown)",
                    data=protocol_report,
                    file_name=f"{timestamp}-KASP_{seq_id}_Protocol.md",
                    mime="text/markdown",
                    use_container_width=True,
                    key="download_protocol_full"
                )
            
            with tab_brief:
                # 生成简化版报告
                brief_report = generate_kasp_protocol_brief(best, seq_id if seq_id else "SNP_Marker")
                st.markdown(brief_report)
                
                # 快速合成单（可直接复制）
                st.markdown("---")
                st.markdown("### 📝 引物合成订单表（可直接复制）")
                
                # 创建可复制的表格
                FAM_TAIL = "gaaggtgaccaagttcatgct"
                HEX_TAIL = "gaaggtcggagtcaacggatt"
                
                order_table = f"""| 引物名称 | 序列 (5'→3') | 长度 | 纯化方式 |
|----------|--------------|------|----------|
| {seq_id}_FAM | {FAM_TAIL}{best['fwd_allele1_core'].upper()} | {len(best['fwd_allele1_full'])} bp | PAGE/HPLC |
| {seq_id}_HEX | {HEX_TAIL}{best['fwd_allele2_core'].upper()} | {len(best['fwd_allele2_full'])} bp | PAGE/HPLC |
| {seq_id}_Common | {best['reverse'].upper()} | {len(best['reverse'])} bp | 脱盐即可 |"""
                
                st.markdown(order_table)
                
                # 纯文本版本方便复制
                st.markdown("#### 📋 纯文本版本（方便复制到合成订单）")
                plain_text = f"""{seq_id}_FAM\t{FAM_TAIL}{best['fwd_allele1_core'].upper()}
{seq_id}_HEX\t{HEX_TAIL}{best['fwd_allele2_core'].upper()}
{seq_id}_Common\t{best['reverse'].upper()}"""
                
                st.code(plain_text, language=None)
                
                st.download_button(
                    label="📥 下载快速参考卡 (Markdown)",
                    data=brief_report,
                    file_name=f"{timestamp}-KASP_{seq_id}_QuickRef.md",
                    mime="text/markdown",
                    use_container_width=True,
                    key="download_protocol_brief"
                )
                
        except ValueError as e:
            st.error(f"❌ 序列解析错误: {e}")
        except Exception as e:
            st.error(f"❌ 设计过程出错: {e}")


def show_regular_pcr_design():
    """常规PCR引物设计页面"""
    st.markdown("### 🧪 常规PCR引物设计")
    
    # 初始化PCR临时保存字段
    if 'temp_pcr_seq_input' not in st.session_state:
        st.session_state['temp_pcr_seq_input'] = None
    if 'temp_pcr_seq_id' not in st.session_state:
        st.session_state['temp_pcr_seq_id'] = None
    if 'temp_pcr_num_pairs' not in st.session_state:
        st.session_state['temp_pcr_num_pairs'] = 5
    if 'temp_pcr_product_min' not in st.session_state:
        st.session_state['temp_pcr_product_min'] = 150
    if 'temp_pcr_product_max' not in st.session_state:
        st.session_state['temp_pcr_product_max'] = 500
    if 'temp_pcr_wheat_mode' not in st.session_state:
        st.session_state['temp_pcr_wheat_mode'] = False
    if 'temp_pcr_avoid_5prime' not in st.session_state:
        st.session_state['temp_pcr_avoid_5prime'] = 40
    if 'temp_pcr_use_target' not in st.session_state:
        st.session_state['temp_pcr_use_target'] = False
    if 'temp_pcr_target_start' not in st.session_state:
        st.session_state['temp_pcr_target_start'] = 50
    if 'temp_pcr_target_end' not in st.session_state:
        st.session_state['temp_pcr_target_end'] = 200
    
    # 显示临时输入恢复提示
    has_temp_pcr_input = st.session_state['temp_pcr_seq_input'] is not None
    if has_temp_pcr_input:
        st.info(f"💫 检测到未保存的编辑内容 | [{st.session_state['temp_pcr_seq_id'] or '未命名'}]")
        temp_col1, temp_col2 = st.columns([1, 1])
        with temp_col1:
            if st.button("✏️ 恢复编辑", use_container_width=True, key="restore_temp_pcr_input"):
                pass
        with temp_col2:
            if st.button("🗑️ 清除草稿", use_container_width=True, key="clear_temp_pcr_input"):
                st.session_state['temp_pcr_seq_input'] = None
                st.session_state['temp_pcr_seq_id'] = None
                st.rerun()
    
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
        # 判断使用临时保存值还是示例
        default_seq = st.session_state['temp_pcr_seq_input'] if st.session_state['temp_pcr_seq_input'] else example_seq
        default_id = st.session_state['temp_pcr_seq_id'] or "My_Gene"
        
        seq_input = st.text_area(
            "输入序列",
            value=default_seq,
            height=200,
            help="输入纯碱基序列(A/T/G/C)",
            key="pcr_seq_input_area"
        )
        # 自动保存到临时缓冲区
        st.session_state['temp_pcr_seq_input'] = seq_input
        
        seq_id = st.text_input("序列ID（可选）", value=default_id, key="regular_seq_id")
        # 自动保存序列ID
        st.session_state['temp_pcr_seq_id'] = seq_id
    
    with col2:
        st.markdown("**参数设置**")
        num_pairs = st.slider("生成引物对数", 3, 10, st.session_state.get('temp_pcr_num_pairs', 5), key="regular_num")
        st.session_state['temp_pcr_num_pairs'] = num_pairs
        
        st.markdown("**产物大小范围**")
        product_min = st.number_input("最小产物(bp)", 100, 500, st.session_state.get('temp_pcr_product_min', 150), key="pcr_min_product")
        st.session_state['temp_pcr_product_min'] = product_min
        product_max = st.number_input("最大产物(bp)", 200, 1000, st.session_state.get('temp_pcr_product_max', 500), key="pcr_max_product")
        st.session_state['temp_pcr_product_max'] = product_max
        
        # 小麦模式
        st.markdown("---")
        wheat_mode = st.checkbox("🌾 小麦特异性模式", value=st.session_state.get('temp_pcr_wheat_mode', False),
                                  help="针对小麦A/B/D同源基因优化，避开5'端保守区", key="pcr_wheat_mode")
        st.session_state['temp_pcr_wheat_mode'] = wheat_mode
        
        if wheat_mode:
            st.info("""**小麦模式已启用：**
- 自动避开5'端保守区域
- 优先在3'端/UTR区设计引物
- 评估同源基因特异性""")
            
            with st.expander("🌾 小麦参数设置"):
                avoid_5prime = st.slider("避开5'端区域(%)", 20, 60, st.session_state.get('temp_pcr_avoid_5prime', 40),
                                         help="避开序列5'端的百分比，该区域同源基因通常高度保守", key="pcr_avoid_5prime")
                st.session_state['temp_pcr_avoid_5prime'] = avoid_5prime
                prefer_3prime = st.checkbox("优先3'端区域", value=True,
                                            help="3'UTR区域通常变异更多，有利于特异性扩增", key="pcr_prefer_3prime")
        else:
            avoid_5prime = 40
            prefer_3prime = True
        
        with st.expander("目标区域（可选）"):
            use_target = st.checkbox("指定目标区域", value=st.session_state.get('temp_pcr_use_target', False), key="pcr_use_target")
            st.session_state['temp_pcr_use_target'] = use_target
            if use_target:
                target_start = st.number_input("起始位置", 1, 10000, st.session_state.get('temp_pcr_target_start', 50), key="pcr_target_start")
                st.session_state['temp_pcr_target_start'] = target_start
                target_end = st.number_input("结束位置", 1, 10000, st.session_state.get('temp_pcr_target_end', 200), key="pcr_target_end")
                st.session_state['temp_pcr_target_end'] = target_end
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
            
            # 生物学验证提示
            st.markdown("---")
            st.markdown("### 🧪 生物学验证（推荐）")
            
            with st.expander("📍 Ensembl Plants BLAST验证步骤", expanded=False):
                st.markdown("""
**目标：** 确认你的引物序列在目标物种上有完美匹配，不在其他物种或非目标区域产生高匹配。

**步骤：**

1. **打开工具**
   - 访问 [Ensembl Plants BLAST](https://plants.ensembl.org/Multi/Tools/Blast?db=core)
   - 或使用 [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/)

2. **输入序列**
   - 复制你的正向或反向引物序列
   - 粘贴到BLAST输入框

3. **关键设置**
   - **Search against:** 选择目标物种及其参考基因组版本
   - **Program:** megablast (对于相同或相似的序列)
   - 保持其他参数默认

4. **点击 Run BLAST**

5. **分析结果**
   
   ✅ **理想情况：**
   - 看到一个 **100% Identity** 的匹配，对应你的 **目标基因组/染色体**
   - 其他可能的匹配度 <90% 或中间有大的 Gap
   - 说明你的引物特异性良好 ✓

   ⚠️ **问题情况：**
   - 在 **非目标区域** 看到 **99%-100% 匹配** 
   - 说明这条序列存在 **同源干扰**，不能用
   - **解决方案：** 修改引物设计或选择其他方案

**小麦特例：** 如果你针对小麦设计，建议用 [PolyMarker](http://polymarker.tgac.ac.uk/) 工具，自动对 A/B/D 三个基因组进行 BLAST。

**注意：** 本工具不直接实现BLAST功能，但强烈建议在使用引物前进行此验证。
                """)
            
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
    
    # 初始化分析页面临时保存字段
    if 'temp_analysis_mode' not in st.session_state:
        st.session_state['temp_analysis_mode'] = 0
    if 'temp_primer_input' not in st.session_state:
        st.session_state['temp_primer_input'] = ""
    if 'temp_primer_name' not in st.session_state:
        st.session_state['temp_primer_name'] = ""
    if 'temp_primer_type' not in st.session_state:
        st.session_state['temp_primer_type'] = "常规PCR"
    if 'temp_check_wheat' not in st.session_state:
        st.session_state['temp_check_wheat'] = False
    if 'temp_primer1_input' not in st.session_state:
        st.session_state['temp_primer1_input'] = ""
    if 'temp_primer2_input' not in st.session_state:
        st.session_state['temp_primer2_input'] = ""
    if 'temp_primer1_name' not in st.session_state:
        st.session_state['temp_primer1_name'] = ""
    if 'temp_primer2_name' not in st.session_state:
        st.session_state['temp_primer2_name'] = ""
    if 'temp_kasp_allele1_input' not in st.session_state:
        st.session_state['temp_kasp_allele1_input'] = ""
    if 'temp_kasp_allele2_input' not in st.session_state:
        st.session_state['temp_kasp_allele2_input'] = ""
    if 'temp_kasp_common_input' not in st.session_state:
        st.session_state['temp_kasp_common_input'] = ""
    
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
    analysis_modes = ["单引物分析", "引物对分析", "小麦KASP引物分析"]
    current_mode_index = st.session_state.get('temp_analysis_mode', 0)
    
    analysis_mode = st.radio(
        "选择分析模式",
        range(len(analysis_modes)),
        format_func=lambda i: analysis_modes[i],
        index=current_mode_index,
        horizontal=True,
        label_visibility="collapsed",
        key="analysis_mode_radio"
    )
    # 保存到临时缓冲区
    st.session_state['temp_analysis_mode'] = analysis_mode
    
    if analysis_modes[analysis_mode] == "单引物分析":
        st.markdown("#### 📝 输入引物序列")
        primer_input = st.text_area(
            "引物序列 (5'→3')",
            value=st.session_state.get('temp_primer_input', ""),
            placeholder="例如: ATGCGATCGATCGATCGATCG\n只输入ATGC碱基，自动过滤其他字符",
            height=100,
            key="analysis_primer_input"
        )
        # 自动保存
        st.session_state['temp_primer_input'] = primer_input
        
        primer_name = st.text_input("引物名称（可选）", value=st.session_state.get('temp_primer_name', "My_Primer"), key="analysis_primer_name")
        # 自动保存
        st.session_state['temp_primer_name'] = primer_name
        
        # 分析类型选择
        col1, col2 = st.columns(2)
        with col1:
            primer_type_options = ["常规PCR", "KASP Allele", "KASP Common"]
            current_type_idx = primer_type_options.index(st.session_state.get('temp_primer_type', "常规PCR"))
            primer_type = st.selectbox("引物类型", range(len(primer_type_options)), format_func=lambda i: primer_type_options[i], 
                                      index=current_type_idx, label_visibility="collapsed", key="analysis_primer_type")
            st.session_state['temp_primer_type'] = primer_type_options[primer_type]
        with col2:
            check_wheat = st.checkbox("小麦特异性检测", value=st.session_state.get('temp_check_wheat', False), key="analysis_check_wheat")
            # 自动保存
            st.session_state['temp_check_wheat'] = check_wheat
        
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
        
        # 检测是否有引物对临时输入
        has_temp_pair_input = st.session_state.get('temp_primer1_input') and st.session_state.get('temp_primer2_input')
        if has_temp_pair_input:
            st.info(f"💫 检测到未保存的编辑内容")
            pair_col1, pair_col2 = st.columns([1, 1])
            with pair_col1:
                if st.button("✏️ 恢复编辑", use_container_width=True, key="restore_temp_pair_input"):
                    pass
            with pair_col2:
                if st.button("🗑️ 清除草稿", use_container_width=True, key="clear_temp_pair_input"):
                    st.session_state['temp_primer1_input'] = ""
                    st.session_state['temp_primer2_input'] = ""
                    st.session_state['temp_primer1_name'] = ""
                    st.session_state['temp_primer2_name'] = ""
                    st.rerun()
        
        col1, col2 = st.columns(2)
        with col1:
            fwd_input = st.text_area("正向引物 (Forward 5'→3')", 
                                     value=st.session_state.get('temp_primer1_input', ""),
                                     placeholder="ATGCGATCGATCGATCG",
                                     height=100,
                                     key="pair_fwd_input")
            # 自动保存到临时缓冲区
            st.session_state['temp_primer1_input'] = fwd_input
            
            fwd_name = st.text_input("正向引物名称", value=st.session_state.get('temp_primer1_name', "Forward_Primer"), key="pair_fwd_name")
            # 自动保存
            st.session_state['temp_primer1_name'] = fwd_name
        
        with col2:
            rev_input = st.text_area("反向引物 (Reverse 5'→3')", 
                                     value=st.session_state.get('temp_primer2_input', ""),
                                     placeholder="CGATCGATCGATCGAT",
                                     height=100,
                                     key="pair_rev_input")
            # 自动保存到临时缓冲区
            st.session_state['temp_primer2_input'] = rev_input
            
            rev_name = st.text_input("反向引物名称", value=st.session_state.get('temp_primer2_name', "Reverse_Primer"), key="pair_rev_name")
            # 自动保存
            st.session_state['temp_primer2_name'] = rev_name
        
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
        
        # 检测是否有KASP临时输入
        has_temp_kasp_input = (st.session_state.get('temp_kasp_allele1_input') or 
                               st.session_state.get('temp_kasp_allele2_input') or 
                               st.session_state.get('temp_kasp_common_input'))
        if has_temp_kasp_input:
            st.info(f"💫 检测到未保存的编辑内容")
            kasp_col1, kasp_col2 = st.columns([1, 1])
            with kasp_col1:
                if st.button("✏️ 恢复编辑", use_container_width=True, key="restore_temp_kasp_input"):
                    pass
            with kasp_col2:
                if st.button("🗑️ 清除草稿", use_container_width=True, key="clear_temp_kasp_input"):
                    st.session_state['temp_kasp_allele1_input'] = ""
                    st.session_state['temp_kasp_allele2_input'] = ""
                    st.session_state['temp_kasp_common_input'] = ""
                    st.rerun()
        
        st.info("""**小麦KASP引物组包括：**
- 2条等位基因特异性引物（带FAM/HEX荧光尾巴）
- 1条通用反向引物（Common Primer）""")
        
        allele1_input = st.text_area("Allele 1 引物（完整，含FAM尾巴）", 
                                     value=st.session_state.get('temp_kasp_allele1_input', ""),
                                     height=80,
                                     key="kasp_allele1_area")
        # 自动保存到临时缓冲区
        st.session_state['temp_kasp_allele1_input'] = allele1_input
        
        allele2_input = st.text_area("Allele 2 引物（完整，含HEX尾巴）", 
                                     value=st.session_state.get('temp_kasp_allele2_input', ""),
                                     height=80,
                                     key="kasp_allele2_area")
        # 自动保存到临时缓冲区
        st.session_state['temp_kasp_allele2_input'] = allele2_input
        
        common_input = st.text_area("Common 反向引物", 
                                    value=st.session_state.get('temp_kasp_common_input', ""),
                                    height=80,
                                    key="kasp_common_area")
        # 自动保存到临时缓冲区
        st.session_state['temp_kasp_common_input'] = common_input
        
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
    
    # 显示Primer3状态
    if PRIMER3_AVAILABLE:
        st.success("✅ Primer3-py库已加载 - 使用专业热力学计算")
    else:
        st.warning("⚠️ Primer3-py库未安装 - 使用内置算法。建议安装: `pip install primer3-py`")
    
    st.markdown("""
    ## 关于本工具 (v6.0 Primer3重构版)
    
    本工具使用 **Primer3-py** 库进行专业的引物设计和热力学分析，参考了 **polyoligo-kasp** 的设计理念。
    
    ### 核心特性
    - 🔬 **Primer3引擎**: 使用SantaLucia最近邻法精确计算Tm值
    - 🧬 **热力学分析**: 专业的发夹、二聚体dG/Tm计算
    - 🌾 **小麦模式**: 针对六倍体小麦的五大忌检测
    - 📊 **综合评分**: 多维度引物质量评估
    
    ---
    
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
    
    用于设计普通PCR扩增引物对。优先使用Primer3引擎，失败时回退到手动算法。
    
    ### 设计原则
    - 引物长度: 18-25 bp
    - Tm值: 55-68°C (SantaLucia法)
    - GC含量: 40-60%
    - 产物大小: 可自定义
    
    ### 热力学参数 (Primer3)
    - 单价阳离子: 50 mM
    - 二价阳离子: 1.5 mM (Mg²⁺)
    - dNTP浓度: 0.6 mM
    - DNA浓度: 250 nM
    
    ### 注意事项
    - 避免引物3'端自身互补 (发夹 Tm < 45°C)
    - 避免引物对之间形成二聚体 (dG > -9 kcal/mol)
    - 两条引物Tm差异应 <2°C
    
    ---
    
    ## 常见问题
    
    **Q: 为什么设计不出引物？**
    A: 可能原因：
    1. 序列太短
    2. GC含量过高或过低
    3. 序列中有过多重复
    4. 没有找到满足热力学条件的候选
    
    **Q: 如何提高引物质量？**
    A: 
    1. 提供更长的侧翼序列
    2. 选择GC含量适中的区域
    3. 调整参数设置
    4. 安装primer3-py库获得更精确的计算
    
    **Q: Primer3和内置算法有什么区别？**
    A: Primer3使用经过验证的SantaLucia热力学参数，Tm计算误差通常<2°C，
       并提供精确的发夹/二聚体dG和Tm值。内置算法是简化版本。
    """)


# ==================== 主程序 ====================

def main():
    # 初始化页面选择状态
    if 'page' not in st.session_state:
        st.session_state['page'] = "🏠 首页"
    
    # 侧边栏导航
    st.sidebar.markdown("## 🧬 引物设计工具")
    st.sidebar.markdown("**v7.0 Web版 (Primer3)**")
    st.sidebar.markdown("---")
    
    # 使用session_state管理页面选择，支持程序内跳转
    page_options = ["🏠 首页", "🔬 KASP引物设计", "🧪 常规PCR引物设计", "🔍 引物分析", "📖 帮助文档"]
    current_index = page_options.index(st.session_state['page']) if st.session_state['page'] in page_options else 0
    
    # Radio选择，不使用key自动binding，而是手动更新
    selected_index = st.sidebar.radio(
        "选择功能",
        range(len(page_options)),
        format_func=lambda i: page_options[i],
        index=current_index,
        label_visibility="collapsed"
    )
    
    # 更新session_state
    st.session_state['page'] = page_options[selected_index]
    page = st.session_state['page']
    
    st.sidebar.markdown("---")
    
    # 显示Primer3状态
    show_primer3_status()
    
    # 显示KASP核心功能说明
    show_kasp_features_info()
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
    <small>
    
    **关于本工具 v7.0**
    
    本工具用于设计KASP基因分型引物和常规PCR引物。
    
    **核心算法: Primer3-py**
    - Tm计算 (SantaLucia法)
    - 发夹结构检测
    - 二聚体风险评估
    - 专业引物设计引擎
    
    **v7.0 新功能:**
    - 🔬 LGC标准人工错配
    - 🆘 AT-rich序列救援模式
    - ⚖️ ASP-Common Tm平衡
    
    </small>
    """, unsafe_allow_html=True)
    
    # 页面路由
    if page == "🏠 首页":
        st.markdown('<p class="main-header">🧬 引物设计工具</p>', unsafe_allow_html=True)
        st.markdown('<p class="sub-header">KASP & 常规PCR 引物设计平台 v6.0 (Primer3-py)</p>', unsafe_allow_html=True)
        
        # 显示Primer3状态卡片
        if PRIMER3_AVAILABLE:
            st.success(f"🔬 **Primer3-py v{PRIMER3_VERSION}** 已加载 - 使用专业热力学算法")
        else:
            st.warning("⚠️ Primer3-py未安装，正在使用内置算法。建议运行 `pip install primer3-py` 获得更精确的计算结果。")
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("""
            ### 🔬 KASP引物设计
            
            针对SNP位点设计KASP基因分型引物
            
            - ✅ 自动添加FAM/HEX荧光尾巴
            - ✅ **LGC标准人工错配** (n-3位置)
            - ✅ **救援模式** (AT-rich序列)
            - ✅ **Tm平衡优化** (ASP-Common匹配)
            - ✅ 多方案评分排序
            - ✅ 🌾 小麦模式 (五大忌检测)
            
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
        
        # v7.0 新功能说明
        st.markdown("---")
        st.markdown("### ✨ v7.0 新功能亮点")
        
        feat_col1, feat_col2, feat_col3, feat_col4 = st.columns(4)
        
        with feat_col1:
            st.info("""
            **🔬 LGC人工错配**
            
            在ASP引物n-3位置引入deliberate mismatch，根据SNP强度选择最佳错配碱基。
            """)
        
        with feat_col2:
            st.info("""
            **🆘 救援模式**
            
            自动检测AT-rich序列，放宽参数确保设计成功，告别空结果。
            """)
        
        with feat_col3:
            st.info("""
            **⚖️ Tm平衡**
            
            Common引物智能匹配ASP的Tm值，确保PCR效率均衡。
            """)
        
        with feat_col4:
            st.info("""
            **🔍 智能诊断**
            
            设计失败时自动分析原因（GC、发夹、复杂度），给出具体建议。
            """)
        
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
