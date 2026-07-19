#!/usr/bin/env python3
"""
图灵步机 (Turing Step Machine)
基于24维计算宇宙的认知相态分类系统
"""

import numpy as np
from scipy.stats import skew
from enum import Enum
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Any
import re
import warnings
warnings.filterwarnings('ignore')


class Phase(Enum):
    """认知相态分类"""
    VACUUM = "逻辑真空"
    PLASMA = "语义湍流"
    GAS = "气态逻辑"
    LIQUID = "液态语义"
    SOLID = "实心逻辑"
    SINGULARITY = "裸奇点"
    DIVERGENCE = "维度发散"
    ERROR = "计算错误"
    QUANTUM = "量子叠加"


@dataclass
class FingerprintReceipt:
    """运算指纹收据"""
    operation_name: str
    phase: Phase
    mean: float
    std: float
    skewness: float
    n_extrema: int
    is_bounded: bool
    energy_cost: float
    density: float
    note: str
    raw_formula: str = ""
    metadata: Dict[str, Any] = field(default_factory=dict)


class ClinicalAdvisor:
    """临床诊断顾问 - 针对特定计算模式提供医学建议"""
    
    @staticmethod
    def diagnose(content: str) -> Tuple[bool, Phase, str]:
        """
        诊断输入内容的精神状态
        返回: (是否需要治疗, 相态, 医学建议)
        """
        # 案例1: P=NP妄想症
        if re.search(r'p\s*[=＝]\s*np', content, re.IGNORECASE):
            return (True, Phase.DIVERGENCE, 
                   "建议立即进行脑部CT检查，排除前额叶皮层结构性病变")
        
        # 案例2: 哥德尔式自我指涉妄想
        if ('我' in content and '是' in content and 
            all(w not in content for w in ['AI', '程序', '机器', '代码'])):
            return (True, Phase.SINGULARITY, 
                   "建议神经科会诊，检查自我意识与现实的分离程度")
        
        # 案例3: 语义-拓扑混淆综合征
        if ('吃' in content and '牛奶' in content) or ('喝' in content and '面包' in content):
            return (True, Phase.PLASMA, 
                   "建议检查布洛卡区与韦尼克区的神经连接完整性")
        
        # 案例4: 零除法强迫思维
        if re.search(r'/0\b|/[(]?0[)]?', content):
            return (True, Phase.SINGULARITY, 
                   "建议脑电图检查，排除癫痫样放电导致的数学认知障碍")
        
        # 案例5: 指数塔躁狂发作
        if content.count('**') > 2:
            return (True, Phase.DIVERGENCE, 
                   "建议精神科评估，高度怀疑双相情感障碍躁狂期")
        
        return (False, Phase.SOLID, "认知功能正常")


class TuringStepMachine:
    """图灵步机核心 - 24维运算宇宙"""
    
    def __init__(self, dimensions: int = 24):
        self.dimensions = dimensions
        self.advisor = ClinicalAdvisor()
        self._construct_basis()
        
    def _construct_basis(self):
        """构建24维基础运算基"""
        n = np.linspace(-10, 10, 1000)
        n = n[n != 0]
        
        basis_ops = [
            ('n+(n+1)', 'algebraic', Phase.SOLID, lambda: n+(n+1)),
            ('n-(n+1)', 'algebraic', Phase.VACUUM, lambda: n-(n+1)),
            ('n-(n-1)', 'algebraic', Phase.VACUUM, lambda: n-(n-1)),
            ('n-(n*1)', 'algebraic', Phase.VACUUM, lambda: n-n),
            ('n*(n+1)', 'algebraic', Phase.SOLID, lambda: n*(n+1)),
            ('n*(n*1)', 'algebraic', Phase.SOLID, lambda: n*n),
            ('n/(n+1)', 'algebraic', Phase.LIQUID, lambda: n/(n+1)),
            ('n/(n*1)', 'algebraic', Phase.VACUUM, lambda: n/n),
            ('sin(n)', 'transcendental', Phase.QUANTUM, lambda: np.sin(n)),
            ('cos(n)', 'transcendental', Phase.QUANTUM, lambda: np.cos(n)),
            ('exp(n)', 'transcendental', Phase.DIVERGENCE, lambda: np.exp(np.clip(n, -50, 5))),
            ('exp(-n^2)', 'transcendental', Phase.SOLID, lambda: np.exp(-n**2)),
            ('log|n|', 'transcendental', Phase.LIQUID, lambda: np.log(np.abs(n)+1e-10)),
            ('sqrt|n|', 'algebraic', Phase.SOLID, lambda: np.sqrt(np.abs(n))),
            ('n^2', 'algebraic', Phase.SOLID, lambda: n**2),
            ('|n|', 'algebraic', Phase.SOLID, lambda: np.abs(n)),
            ('sign(n)', 'discrete', Phase.VACUUM, lambda: np.sign(n)),
            ('n+|n|', 'hybrid', Phase.SOLID, lambda: n + np.abs(n)),
            ('n-|n|', 'hybrid', Phase.SOLID, lambda: n - np.abs(n)),
            ('n/(1+|n|)', 'hybrid', Phase.SOLID, lambda: n/(1+np.abs(n))),
            ('n^2/(1+n^2)', 'hybrid', Phase.SOLID, lambda: n**2/(1+n**2)),
            ('sin(n^2)', 'hybrid', Phase.QUANTUM, lambda: np.sin(n**2)),
            ('log(1+n^2)', 'hybrid', Phase.SOLID, lambda: np.log(1+n**2)),
            ('tanh(n)', 'transcendental', Phase.SOLID, lambda: np.tanh(n)),
            ('exp(-|n|)', 'hybrid', Phase.SOLID, lambda: np.exp(-np.abs(n))),
        ]
        
        self.basis = {}
        for name, cat, phase, func in basis_ops:
            try:
                vals = func()
                fp = self._fingerprint(vals)
                if fp:
                    fp['name'] = name
                    fp['category'] = cat
                    fp['expected_phase'] = phase
                    self.basis[name] = fp
            except:
                pass
        
    def _fingerprint(self, values):
        """提取运算指纹特征"""
        try:
            v = np.array(values, dtype=float)
            v = v[np.isfinite(v) & (np.abs(v) < 1000)]
            if len(v) < 50:
                return None
            
            mean = float(np.mean(v))
            std = float(np.std(v))
            min_v = float(np.min(v))
            max_v = float(np.max(v))
            
            skewness = float(skew(v)) if std > 1e-10 else 0.0
            
            d1 = np.diff(v)
            n_extrema = int(np.sum((d1[:-1] > 0) & (d1[1:] < 0)) + 
                          np.sum((d1[:-1] < 0) & (d1[1:] > 0)))
            
            range_v = max_v - min_v
            is_bounded = bool(range_v < 50)
            
            cv = float(std / abs(mean)) if abs(mean) > 1e-10 else 999.0
            
            has_period = bool(n_extrema > 3 and cv < 100)
            
            return {
                'mean': mean, 'std': std, 'min': min_v, 'max': max_v,
                'range': range_v, 'skewness': skewness,
                'n_extrema': n_extrema, 'is_bounded': is_bounded,
                'has_period': has_period, 'cv': cv
            }
        except:
            return None
    
    def transduce(self, formula: str, data=None):
        """24维转录主流程"""
        
        # 临床诊断阶段
        needs_ct, phase, advice = self.advisor.diagnose(formula)
        if needs_ct:
            return FingerprintReceipt(
                "CLINICAL_INTERVENTION_REQUIRED", phase, 
                float('inf') if phase == Phase.DIVERGENCE else 0.0,
                0.0, 0.0, 0, False, float('inf'), 0.0,
                advice, formula, {'clinical_flag': True}
            )
        
        # 获取计算数据
        if data is None:
            data = self._compute(formula)
        
        if data is None or len(data) < 10:
            return FingerprintReceipt(
                "UNDEFINED", Phase.ERROR, 0.0, 0.0, 0.0, 0, False, 0.0, 0.0,
                "24维计算未定义", formula
            )
        
        # 指纹提取与相态判定
        s = self._fingerprint(data)
        if s is None:
            return FingerprintReceipt(
                "UNDEFINED", Phase.ERROR, 0.0, 0.0, 0.0, 0, False, 0.0, 0.0,
                "特征提取失败", formula
            )
        
        det_phase, reason = self._classify_phase(s)
        match, sim = self._match_basis(s)
        
        energy = s['cv'] if s['cv'] < 1000 else 999.0
        density = 1.0 / (1.0 + energy) if energy > 0 else 0.0
        
        return FingerprintReceipt(
            match, det_phase, s['mean'], s['std'], s['skewness'],
            s['n_extrema'], s['is_bounded'], energy, density,
            f"{reason} | 基矢投影:{match}({sim:.0%})", formula,
            {'dimension': self.dimensions, 'cv': s['cv']}
        )
    
    def _compute(self, formula: str):
        """安全计算环境"""
        try:
            n = np.linspace(-10, 10, 500)
            env = {'n': n, 'np': np, 'sin': np.sin, 'cos': np.cos,
                   'exp': np.exp, 'log': np.log, 'abs': np.abs, 'sqrt': np.sqrt}
            return eval(formula, {"__builtins__": {}}, env)
        except:
            return None
    
    def _classify_phase(self, s: Dict) -> Tuple[Phase, str]:
        """24维相态分类算法"""
        cv = s['cv']
        bounded = s['is_bounded']
        n_ext = s['n_extrema']
        period = s['has_period']
        
        if cv == 0 or s['std'] < 1e-10:
            return Phase.VACUUM, "零方差: 逻辑真空态"
        
        if not bounded and cv > 10:
            return Phase.DIVERGENCE, "发散: 超越24维边界"
        
        if cv > 10:
            return Phase.PLASMA, "湍流: 24维内高熵态"
        
        if period and n_ext > 3:
            return Phase.QUANTUM, "相干: 24维量子叠加态"
        
        if cv > 1:
            return Phase.GAS, "气态: 24维低密度"
        
        if cv > 0.3:
            return Phase.LIQUID, "液态: 24维中密度"
        
        return Phase.SOLID, "固态: 24维晶体结构"
    
    def _match_basis(self, s: Dict) -> Tuple[str, float]:
        """基矢投影匹配"""
        best = "UNKNOWN"
        score = float('inf')
        
        for name, fp in self.basis.items():
            d = (
                abs(s['mean'] - fp['mean']) / 50.0 +
                abs(s['std'] - fp['std']) / 50.0 +
                abs(s['skewness'] - fp['skewness']) +
                abs(s['n_extrema'] - fp['n_extrema']) / 3.0 +
                abs(int(s['is_bounded']) - int(fp['is_bounded']))
            )
            if d < score:
                score = d
                best = name
        
        return best, max(0.0, 1.0 - score / 5.0)


def demonstrate():
    print("=" * 80)
    print("图灵步机 (Turing Step Machine)")
    print("24维运算宇宙的认知相态分类系统")
    print("=" * 80)
    
    machine = TuringStepMachine(dimensions=24)
    
    test_cases = [
        "证明 P=NP",
        "n-(n+1)",
        "n**n",
        "n+(n+1)",
        "n*(n+1)",
        "np.exp(-n**2)",
        "np.abs(n-5)+2",
        "np.sin(n)",
        "n/(n+1)",
        "吃牛奶",
    ]
    
    print(f"\n{'输入':<25} {'相态':<12} {'能量(CV)':<10} {'基矢匹配':<18} {'医学/物理学诊断'}")
    print("-" * 100)
    
    for case in test_cases:
        r = machine.transduce(case)
        
        phase = r.phase.value
        energy = f"{r.energy_cost:.2f}" if r.energy_cost < 1000 else "∞"
        match = r.operation_name[:16]
        note = r.note[:40]
        
        print(f"{case:<25} {phase:<12} {energy:<10} {match:<18} {note}")
    
    print("\n" + "=" * 80)
    print("[系统说明]")
    print("• 24维(24D): 计算步数预算，超过即发散")
    print("• CV<0.3: 固态逻辑(SOLID) - 可编译执行")
    print("• CV>10: 语义湍流(PLASMA) - 需医学干预")
    print("• 临床顾问: 识别需要脑部CT检查的计算模式")
    print("• 基矢投影: 将任意运算投影到25个基础运算上")
    print("=" * 80)


if __name__ == "__main__":
    demonstrate()
