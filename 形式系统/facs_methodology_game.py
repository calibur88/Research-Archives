# ============================================================
# FACS × 方法博弈论 互锁形式化演示
# Formalized Adjudication Calculus System v1.3
# ×
# 元方法博弈论 v3.7
# ============================================================

from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Set, Dict, List, Optional, Callable, Tuple
from collections import deque
import hashlib
import json

# -----------------------------------------------------------
# 一、FACS 核心：基态与常量
# -----------------------------------------------------------

class State(Enum):
    """L1 八种基态"""
    NABLA = "∇"      # 探索态 (L2)
    DIAMOND = "◇"    # 待验态
    SQUARE = "□"     # 临时确认态
    BLACK_DIAMOND = "◆"  # 永久锚定态
    CIRCLE = "○"     # 修订锁定态
    BULLET = "●"     # 回滚痕
    BOT = "⊥"        # 证伪坟场
    ARROW_CIRCLE = "↻"  # 待回提请

class L0Constant(Enum):
    """L0 超验基底常量（一次性注入，全局只读）"""
    PHI = "Φ"        # 语法生成器
    XI = "Ξ"         # 依赖拓扑公理（DAG强制）
    GAMMA = "Γ"      # 自扩展元规则
    PSI = "Ψ"        # 黑箱裁判接口

@dataclass(frozen=True)
class K0_Axiom:  pass  # 排中律在L1内有效
@dataclass(frozen=True)
class K1_Axiom:  pass  # 矛盾律是L1唯一绝对禁止条件
@dataclass(frozen=True) 
class K2_Axiom:
    max_steps: int = 1000  # 证明链长度上限（硬件锚定）

# -----------------------------------------------------------
# 二、命题与依赖图
# -----------------------------------------------------------

@dataclass
class Proposition:
    pid: str                          # 唯一标识
    content: str                      # 命题内容
    state: State = State.DIAMOND      # 默认待验态
    dependencies: Set[str] = field(default_factory=set)  # Dep(C)
    layer: str = "L1"                 # 当前所在层
    history: List[Tuple[str, State]] = field(default_factory=list)  # 状态历史

    def __hash__(self):
        return hash(self.pid)

    def transition(self, new_state: State, reason: str = ""):
        self.history.append((reason, self.state))
        self.state = new_state

class DependencyGraph:
    """有向无环图（DAG）管理器 —— L0 Ξ 的实现"""

    def __init__(self):
        self.nodes: Dict[str, Set[str]] = {}  # pid -> {dep_pids}
        self._cache_valid = True
        self._cycle_free_cache = True

    def add_edge(self, from_pid: str, to_pid: str) -> bool:
        """添加依赖边，若形成环则返回 False"""
        if from_pid not in self.nodes:
            self.nodes[from_pid] = set()
        if to_pid in self.nodes.get(from_pid, set()):
            return True  # 已存在

        # 临时添加并检测环
        self.nodes[from_pid].add(to_pid)
        has_cycle = self._detect_cycle()

        if has_cycle:
            self.nodes[from_pid].remove(to_pid)  # 回滚
            return False
        return True

    def _detect_cycle(self) -> bool:
        """DFS 环检测"""
        WHITE, GRAY, BLACK = 0, 1, 2
        color = {n: WHITE for n in self.nodes}

        # 包含所有被依赖的节点
        for node in list(self.nodes.keys()):
            for dep in self.nodes[node]:
                if dep not in color:
                    color[dep] = WHITE

        def dfs(node):
            color[node] = GRAY
            for neighbor in self.nodes.get(node, set()):
                if color.get(neighbor, BLACK) == GRAY:
                    return True  # 发现回边
                if color.get(neighbor, BLACK) == WHITE and dfs(neighbor):
                    return True
            color[node] = BLACK
            return False

        for node in color:
            if color[node] == WHITE:
                if dfs(node):
                    return True
        return False

    def get_closure(self, pid: str) -> Set[str]:
        """获取依赖闭包"""
        visited = set()
        queue = deque([pid])
        while queue:
            current = queue.popleft()
            if current in visited:
                continue
            visited.add(current)
            for dep in self.nodes.get(current, set()):
                queue.append(dep)
        visited.remove(pid)
        return visited

# -----------------------------------------------------------
# 三、FACS L1 算子系统
# -----------------------------------------------------------

class FACS:
    """
    形式化可检验演算系统
    L0: 只读常量库
    L1: 可检验演算域
    L2: 前理论探索区
    """

    def __init__(self, k2_limit: int = 100):
        # L0: 一次性锚定为 ◆ 的常量
        self.L0_axioms = {
            "K0": (K0_Axiom(), State.BLACK_DIAMOND),
            "K1": (K1_Axiom(), State.BLACK_DIAMOND),
            "K2": (K2_Axiom(max_steps=k2_limit), State.BLACK_DIAMOND),
        }
        self.L0_constants = {c.value: State.BLACK_DIAMOND for c in L0Constant}

        # L1/L2 存储
        self.propositions: Dict[str, Proposition] = {}
        self.dag = DependencyGraph()
        self.L2_pool: Dict[str, Proposition] = {}  # 探索态

        # 黑箱裁判接口 Ψ 的实现
        self.psi_reports: List[dict] = []

    # --- L2 操作 ---
    def inject_L2(self, pid: str, content: str) -> Proposition:
        """注入 L2: 命题初始为 ∇ 态"""
        p = Proposition(pid=pid, content=content, state=State.NABLA, layer="L2")
        self.L2_pool[pid] = p
        return p

    def structural_review(self, pid: str) -> Tuple[bool, str]:
        """
        Ψ.StructuralReview: L2 → L1 晋升审查
        1. 句法合法性
        2. 静态定义域
        3. 操作路径可定位性
        """
        if pid not in self.L2_pool:
            return False, "PID not in L2"
        p = self.L2_pool[pid]

        # 检查1: 句法合法性（简化：内容非空且可哈希）
        if not p.content or len(p.content) < 2:
            return False, "SYNTAX_FAIL: content too short"

        # 检查2: 静态定义域（简化：内容包含明确的 scope 声明）
        if "SCOPE:" not in p.content and "适用域" not in p.content:
            return False, "SCOPE_UNDEFINED: static domain missing"

        # 检查3: 操作路径可定位性（简化：包含可操作路径描述）
        if "PATH:" not in p.content and "操作路径" not in p.content:
            return False, "PATH_UNLOCATABLE: no operational trace"

        return True, "PASS"

    def promote(self, pid: str) -> bool:
        """L2 → L1 晋升"""
        ok, msg = self.structural_review(pid)
        if not ok:
            print(f"  [晋升拒绝] {pid}: {msg}")
            return False
        p = self.L2_pool.pop(pid)
        p.layer = "L1"
        p.state = State.DIAMOND
        self.propositions[pid] = p
        print(f"  [晋升成功] {pid}: ∇ → ◇")
        return True

    # --- L1 算子 ---
    def confirm(self, pid: str) -> bool:
        """确认 □: ◇C → □C"""
        if pid not in self.propositions:
            return False
        p = self.propositions[pid]
        if p.state != State.DIAMOND:
            return False

        # 检查依赖无环（Ξ）
        for dep in p.dependencies:
            if not self.dag.add_edge(pid, dep):
                print(f"  [确认失败] {pid}: CIRCULAR_DEP detected")
                p.transition(State.BOT, "CIRCULAR_DEP")
                return False

        # Ψ.Consistency 检查
        if not self._check_consistency(pid):
            return False

        p.transition(State.SQUARE, "confirm")
        print(f"  [确认] {pid}: ◇ → □")
        return True

    def finalize(self, pid: str) -> bool:
        """永久锚定 ◆: □C → ◆C，写入 L0 只读库"""
        if pid not in self.propositions:
            return False
        p = self.propositions[pid]
        if p.state != State.SQUARE:
            return False

        # 检查：全部依赖为 ◆ 或初始公设
        for dep_id in p.dependencies:
            dep = self.propositions.get(dep_id)
            if dep and dep.state != State.BLACK_DIAMOND:
                if dep_id not in self.L0_axioms:
                    print(f"  [锚定失败] {pid}: 依赖 {dep_id} 非 ◆")
                    return False

        p.transition(State.BLACK_DIAMOND, "finalize")
        # 从动态 DAG 物理脱离，写入 L0 逻辑库
        print(f"  [锚定] {pid}: □ → ◆ (L0)")
        return True

    def retract(self, pid: str) -> bool:
        """级联回滚 ●: □C → ●C，递归撤销依赖闭包，遇 ◆ 抛 AnchorGuard"""
        if pid not in self.propositions:
            return False
        p = self.propositions[pid]
        if p.state not in (State.SQUARE, State.DIAMOND):
            return False

        closure = self.dag.get_closure(pid)
        for cid in closure:
            cp = self.propositions.get(cid)
            if cp:
                if cp.state == State.BLACK_DIAMOND:
                    print(f"  [AnchorGuard] 遇 ◆ {cid}，中断回滚")
                    return False
                cp.transition(State.BULLET, f"cascade_retract_from_{pid}")
                print(f"  [回滚] {cid}: → ●")

        p.transition(State.BULLET, "retract_root")
        print(f"  [回滚根] {pid}: → ●")
        return True

    def falsify(self, pid: str) -> bool:
        """自我证伪 ¬◇: 推出 C → ⊥"""
        if pid not in self.propositions:
            return False
        p = self.propositions[pid]

        if p.state == State.BLACK_DIAMOND:
            print(f"  [静默拒绝] {pid}: ◆ 免疫证伪")
            return False

        if p.state == State.SQUARE:
            self.retract(pid)  # 先强制回滚

        p.transition(State.BOT, "falsify")
        print(f"  [证伪] {pid}: → ⊥")
        return True

    def _check_consistency(self, pid: str) -> bool:
        """Ψ.Consistency: 检查与已确认集无显性冲突"""
        # 简化实现：检查内容是否同时包含 A 和 ¬A
        p = self.propositions[pid]
        content = p.content
        # 模拟：若内容包含 "CONTRADICT" 则判定不一致
        if "CONTRADICT" in content:
            print(f"  [一致性失败] {pid}: explicit contradiction found")
            return False
        return True

    # --- 核心：全域一致性对撞协议 ---
    def collision_protocol(self, s1_pid: str, s2_pid: str, k2_steps: int = 50) -> dict:
        """
        第六节：全域一致性对撞协议
        检查 Th(S₁) ∩ Th(S₂) = ∅ ? 在 K2 步内
        """
        print(f"\n{'='*60}")
        print(f"[对撞协议启动] {s1_pid}  vs  {s2_pid}")
        print(f"{'='*60}")

        s1 = self.propositions.get(s1_pid)
        s2 = self.propositions.get(s2_pid)

        if not s1 or not s2:
            return {"verdict": "ERROR", "reason": "system_not_found"}

        # 阶段1: 句法审查（已晋升到 L1 时完成）
        # 阶段2: 独立内部一致性（沙箱推演）
        print("\n[阶段2] 独立沙箱推演...")

        # 模拟：解析两个系统的"定理"
        th1 = self._extract_theorems(s1)
        th2 = self._extract_theorems(s2)

        # 阶段3: 合并对撞
        print("[阶段3] 合并图构建与推演...")
        intersection = th1 & th2

        # 阶段4: 终局裁决
        report = {
            "s1": s1_pid, "s2": s2_pid,
            "th1_size": len(th1), "th2_size": len(th2),
            "intersection": list(intersection),
            "steps_used": min(k2_steps, len(th1) + len(th2)),
        }

        if intersection:
            # 情况 A: 矛盾
            print(f"[阶段4] 情况 A: 捕获矛盾命题 → 修订锁定 ○")
            s1.transition(State.CIRCLE, "collision_contradiction")
            s2.transition(State.CIRCLE, "collision_contradiction")
            report["verdict"] = "CONTRADICTION"
            report["action"] = "REVISION_REQUIRED"
            # 生成 ◆META
            meta = self._finalize_meta(report)
            report["meta_id"] = meta
        else:
            # 情况 B: 正交
            print(f"[阶段4] 情况 B: 逻辑正交，解除锁定 → □")
            if s1.state == State.CIRCLE:
                s1.transition(State.SQUARE, "collision_orthogonal")
            if s2.state == State.CIRCLE:
                s2.transition(State.SQUARE, "collision_orthogonal")
            report["verdict"] = "ORTHOGONAL"
            report["action"] = "INDEPENDENT_EVOLUTION"
            meta = self._finalize_meta(report)
            report["meta_id"] = meta

        self.psi_reports.append(report)
        return report

    def _extract_theorems(self, p: Proposition) -> Set[str]:
        """模拟句法后承提取"""
        # 从内容中提取 "THEOREM:xxx" 格式的命题
        theorems = set()
        for line in p.content.split("\n"):
            if line.strip().startswith("THEOREM:"):
                theorems.add(line.strip().replace("THEOREM:", "").strip())
        return theorems

    def _finalize_meta(self, report: dict) -> str:
        """Ψ.Finalize: 裁决报告锚定为 ◆META"""
        meta_id = f"META_{hashlib.sha256(json.dumps(report, sort_keys=True).encode()).hexdigest()[:8]}"
        meta_p = Proposition(
            pid=meta_id,
            content=json.dumps(report, ensure_ascii=False),
            state=State.BLACK_DIAMOND,
            layer="L0"
        )
        self.propositions[meta_id] = meta_p
        print(f"  [元锚定] {meta_id}: → ◆META (L0 只读)")
        return meta_id

    def status(self):
        """状态查询"""
        print(f"\n{'─'*50}")
        print("[FACS 全局状态]")
        print(f"L0 常量: {list(self.L0_constants.keys())}")
        print(f"L2 探索区: {list(self.L2_pool.keys())}")
        print(f"L1 命题:")
        for pid, p in self.propositions.items():
            if not pid.startswith("META"):
                print(f"  {pid}: [{p.state.value}] (层:{p.layer})")
        print(f"L0 锚定元数据: {[k for k in self.propositions if k.startswith('META')]}")
        print(f"{'─'*50}")


# -----------------------------------------------------------
# 四、方法博弈论：三层判据体系
# -----------------------------------------------------------

class Layer(Enum):
    """博弈论三层"""
    EXISTENTIAL = "存在层"      # 第一层
    VERIFICATORY = "真值层"     # 第二层
    INVESTIGATIVE = "研究层"    # 第三层

@dataclass
class Theory:
    """理论对象"""
    tid: str
    content: str           # 完整文本
    axioms: List[str] = field(default_factory=list)
    scope: str = ""        # 静态定义域
    path_desc: str = ""    # 操作路径描述
    layer: Layer = Layer.INVESTIGATIVE
    existence_granted: bool = False  # 是否跨过存在性门槛

class MethodologyGame:
    """
    元方法博弈论
    核心规则：存在权与真值权严格分流
    """

    def __init__(self):
        self.theories: Dict[str, Theory] = {}
        self.existence_registry: Set[str] = set()  # 已授予存在权的理论
        self.verdict_history: List[dict] = []

    def submit(self, tid: str, content: str, axioms: List[str], scope: str, path_desc: str):
        """提交新理论到研究层"""
        t = Theory(tid=tid, content=content, axioms=axioms, 
                   scope=scope, path_desc=path_desc, layer=Layer.INVESTIGATIVE)
        self.theories[tid] = t
        print(f"[博弈论] 理论 '{tid}' 提交至研究层")
        return t

    # --- 第一层：存在性门槛审查 ---
    def existential_review(self, tid: str) -> Tuple[bool, List[str]]:
        """
        4.1 存在层审查：四项合取
        (1) 语法一致性
        (2) 可客观翻译性
        (3) 静态定义域归属
        (4) 操作路径可定位性
        """
        if tid not in self.theories:
            return False, ["THEORY_NOT_FOUND"]
        t = self.theories[tid]
        checks = []

        # (1) 语法一致性：检查显性矛盾 A ∧ ¬A
        has_explicit_contra = self._check_explicit_contradiction(t.content)
        if has_explicit_contra:
            checks.append("FAIL_SYNTAX: Explicit contradiction detected")
        else:
            checks.append("PASS_SYNTAX")

        # (2) 可客观翻译性：核心断言可映射到其他符号系统
        translatable = self._check_translatability(t)
        if not translatable:
            checks.append("FAIL_TRANSLATION: Core assertions not mappable")
        else:
            checks.append("PASS_TRANSLATION")

        # (3) 静态定义域归属
        scope_static = len(t.scope) > 0 and "滑移" not in t.scope
        if not scope_static:
            checks.append("FAIL_SCOPE: Static domain undefined or slippery")
        else:
            checks.append("PASS_SCOPE")

        # (4) 操作路径可定位性（有限资源）
        path_ok = len(t.path_desc) > 0 and "无限" not in t.path_desc
        if not path_ok:
            checks.append("FAIL_PATH: Requires infinite resources or unlocatable")
        else:
            checks.append("PASS_PATH")

        all_pass = all(c.startswith("PASS") for c in checks)
        return all_pass, checks

    def _check_explicit_contradiction(self, content: str) -> bool:
        """检查显性矛盾：同一行同时出现 A 和 ¬A"""
        lines = content.split("\n")
        for line in lines:
            pass
        return "EXPLICIT_CONTRADICT" in content

    def _check_translatability(self, t: Theory) -> bool:
        """可客观翻译性：检查是否有明确的映射规则"""
        return "TRANSLATE:" in t.content or "映射规则" in t.content

    def grant_existence(self, tid: str) -> bool:
        """授予存在权：跨域晋升"""
        ok, checks = self.existential_review(tid)
        t = self.theories[tid]

        print(f"\n[博弈论·存在层审查] 理论 '{tid}':")
        for c in checks:
            print(f"  {c}")

        if ok:
            t.layer = Layer.EXISTENTIAL
            t.existence_granted = True
            self.existence_registry.add(tid)
            print(f"  → 存在权授予: '{tid}' 跨越存在性门槛")
            print(f"  → 本体论不可追溯权生效")
            return True
        else:
            print(f"  → 存在权拒绝: '{tid}' 停留在研究层")
            return False

    # --- 第二层：真值层审查（由外部系统如 FACS 执行）---
    def verificatory_note(self, tid: str, verdict: str, details: str):
        """记录真值层裁决（不追溯撤销存在权）"""
        record = {
            "tid": tid,
            "layer": "VERIFICATORY",
            "verdict": verdict,
            "details": details,
            "note": "存在权不可追溯撤销"
        }
        self.verdict_history.append(record)
        print(f"[博弈论·真值层记录] '{tid}': {verdict} — 存在权保持完整")

    def report(self):
        print(f"\n{'='*60}")
        print("[方法博弈论·全局报告]")
        print(f"研究层理论: {[t for t, v in self.theories.items() if v.layer == Layer.INVESTIGATIVE]}")
        print(f"存在层理论: {list(self.existence_registry)}")
        print(f"真值层裁决记录: {len(self.verdict_history)} 条")
        for v in self.verdict_history:
            print(f"  {v['tid']}: {v['verdict']}")
        print(f"{'='*60}")


# -----------------------------------------------------------
# 五、互锁演示主函数
# -----------------------------------------------------------

def demo_mutual_lock():
    """
    演示 FACS 与 方法博弈论 的互锁闭环
    """
    print("╔" + "═"*78 + "╗")
    print("║" + " FACS × 方法博弈论 · 互锁形式化演示 ".center(76) + "║")
    print("║" + " Mutual Lock: Existence ⟷ Consistency ".center(76) + "║")
    print("╚" + "═"*78 + "╝")

    facs = FACS(k2_limit=100)
    mg = MethodologyGame()

    # 阶段 B: 博弈论审查 FACS 的存在权
    print("\n" + "▓"*80)
    print("▓ 阶段 B: 博弈论审查 FACS 的存在权")
    print("▓"*80)

    facs_content = """
SYSTEM: FACS_v1.3
TYPE: 自举判定系统
SCOPE: 适用于一切具备可定义依赖关系的形式化命题系统
TRANSLATE: 全部术语指涉定义严格锁定于本文框架内部
PATH: 任意命题在 K2 步内完成状态转移穷举
AXIOMS: Φ Ξ Γ Ψ K0 K1 K2
"""
    mg.submit("FACS_v1.3", facs_content, ["Φ","Ξ","Γ","Ψ","K0","K1","K2"],
              "形式化命题系统", "有限步图算法")
    mg.grant_existence("FACS_v1.3")

    # 阶段 C: FACS 构造博弈论系统对
    print("\n" + "▓"*80)
    print("▓ 阶段 C: FACS 构造博弈论及其逆系统")
    print("▓"*80)

    mg_content = """
SYSTEM: MethodologyGame
SCOPE: 元方法论
PATH: 有限资源审查
THEOREM: existence_precedes_verification
THEOREM: ontological_non_retroactivity
"""
    not_mg_content = """
SYSTEM: NOT_MethodologyGame
SCOPE: 元方法论逆
PATH: 有限资源审查
THEOREM: verification_precedes_existence
THEOREM: ontological_retroactivity_enabled
"""
    facs.inject_L2("MG_S1", mg_content)
    facs.inject_L2("MG_S2", not_mg_content)
    facs.promote("MG_S1")
    facs.promote("MG_S2")
    facs.confirm("MG_S1")
    facs.confirm("MG_S2")
    facs.status()

    # 阶段 D: 对撞协议
    print("\n" + "▓"*80)
    print("▓ 阶段 D: 全域一致性对撞协议")
    print("▓"*80)
    report = facs.collision_protocol("MG_S1", "MG_S2", k2_steps=100)

    # 阶段 E: 裁决回写
    print("\n" + "▓"*80)
    print("▓ 阶段 E: 裁决回写与方法论闭环")
    print("▓"*80)

    if report["verdict"] == "CONTRADICTION":
        print("\n[结果] 两系统闭包相交 —— 矛盾 detected")
        mg.verificatory_note("MethodologyGame", "COLLISION_CONTRADICTION", "需修订")
    else:
        print("\n[结果] 两系统逻辑正交 —— 无交集")
        mg.verificatory_note("MethodologyGame", "COLLISION_ORTHOGONAL", "一致性通过")

    # 阶段 F: 元数据锚定
    print("\n" + "▓"*80)
    print("▓ 阶段 F: 元数据锚定 (◆META)")
    print("▓"*80)
    print(f"\n  裁决报告已锚定为 L0 只读: {report.get('meta_id', 'N/A')}")

    # 互锁结构总结
    print("\n" + "╔" + "═"*78 + "╗")
    print("║" + " 互锁闭环验证 ".center(76) + "║")
    print("╠" + "═"*78 + "╣")
    print("║  ① 博弈论授予 FACS 存在权 (第一层)".ljust(76) + "║")
    print("║     ↓ 存在权不可追溯，不依赖 FACS 的一致性结果".ljust(76) + "║")
    print("║  ② FACS 审查博弈论一致性 (第二层)".ljust(76) + "║")
    print("║     ↓ 一致性审查不撤销博弈论的存在权".ljust(76) + "║")
    print("║  ③ 双向互锁，但 DAG 中无环".ljust(76) + "║")
    print("║  ④ 终局裁决: ◆META 锚定，不可上诉".ljust(76) + "║")
    print("╚" + "═"*78 + "╝")

    facs.status()
    mg.report()
    return facs, mg, report


if __name__ == "__main__":
    demo_mutual_lock()
