# NBG 内禀自守量子公理系统

## 量子性的纯代数起源：无需外部参数

**作者**：CH.HY  
**所属机构**：独立研究者  
**邮箱**：jiuxin303@qq.com  
**ORCID**：0009-0003-6134-3736  
**日期**：2026-07-15

---

## 摘要

我们在冯·诺依曼–伯奈斯–哥德尔（NBG）类论的框架内，给出了量子理论的一个完备公理化基础。量子性并非通过固定的普朗克常数或测量假设来强加，而是源于作用在自守形式上的类傅里叶变换在所有序数层级上递归刚性的破缺。柯西–施瓦茨等式被等同于全局递归相容性，这迫使基本变量的对易子消失，谱为纯点状（紧致）。在任何序数层级上，递归比例性的单次失效都会立即触发严格不等式、非零对易子、连续谱以及非紧性。代数普朗克常数被几何地等同为自守线丛的跑动第一陈类，分解为一个裸整数拓扑不变量和一个与尺度相关的真空涨落修正，直接对应 QCD 拓扑磁化率和手征对称破缺序参量。该系统不容许中间态；经典/量子二分法由序数递归律逻辑推出。所有传统欧氏 QFT 公理均为本框架在经典不动点处的特例极限。本文进一步证明，氢原子的离散–连续谱结构是序数相变的物理原型，斯特恩–盖拉赫实验是投影几何不可交换性的判决性演示，波粒二象性是表象相变的两面，而迹反常则是逆梅林变换在有限截断下的代数残余。

---

## 零、基础定义（唯一底层语言）

### 0.1 空间基底

取正实数乘法群 $\mathbb{R}_+^\times$，配备 Haar 不变测度：

$$
d^\times x = \frac{dx}{x}
$$

平方可积空间：

$$
\mathcal{H} = L^2(\mathbb{R}_+^\times, d^\times x)
$$

对数酉变换 $X = \log x$ 给出空间酉等价：

$$
\mathcal{H} \cong L^2(\mathbb{R}, dX)
$$

### 0.2 核心规范算子

定义标度对称生成算子：

$$
A = -(xp + px)
$$

利用正则恒等式 $px = xp - i$，算子恒等式严格成立：

$$
xp + px = 2xp - i
$$

### 0.3 梅林表象（本体系标准表象）

定义梅林变换：

$$
\mathcal{M}[f](s) = \int_0^\infty x^{s-1} f(x)\, dx
$$

对基本算子有固定像：

- $\mathcal{M}[xp](s) = is$
- $\mathcal{M}[px](s) = i(s-1)$

二者在梅林空间均为纯乘法算子。

**关键恒等式（本体系根基）**：

$$
[px,\; xp] = [xp - i,\; xp] = 0
$$

**结论**：梅林表象下位置、动量复合算子天然对易。时域非对易是逆梅林变换投影残留，不是算子内禀属性。

---

## 一、公理 0：基空间与对象

设 $\Gamma$ 为 NBG 中的一个类离散群，$X$ 为一个类对称空间。量子态空间定义为广义自守空间：

$$
\mathcal{A}!ut(\Gamma \backslash X, k)
$$

内禀变量类 $\mathfrak{X}, \mathfrak{Y}$ 是权为 $k$ 的自守形式，满足：

$$
\mathfrak{X}(\gamma z) = j(\gamma, z)^k \mathfrak{X}(z), \quad \forall \gamma \in \Gamma.
$$

它们生成一个算子环 $\mathcal{A}$。类上的基本双线性内积由类求和 $\langle \cdot, \cdot \rangle := \sum$ 定义，并诱导出范数 $\|\cdot\|$。

*注记*：类求和在自守类上是良定义的。在情形 A（紧算子，纯点谱）中，该和是可数收敛的；在情形 B 中，我们引入尺度加权内积 $\langle \cdot, \cdot \rangle_\mu$ 来正则化长程发散类求和。

---

## 二、公理 I：NBG 递归柯西刚性判据（基本定律）

### 类空间上的柯西不等式

$$
\|\mathfrak{X}\|^2 \|\mathfrak{Y}\|^2 \ge |\langle \mathfrak{X}, \mathfrak{Y} \rangle|^2
$$

### 递归等式的定义

设 $\mathcal{F}$ 为类傅里叶变换，$\mathcal{F}^n$ 为序数层级 $n$ 上的第 $n$ 次递归傅里叶映射。对于极限序数 $\lambda = \sup_{\alpha<\lambda} \alpha$，定义投射极限 $\mathcal{F}^\lambda := \varprojlim_{\alpha<\lambda} \mathcal{F}^\alpha$，并假定该极限在类 $\mathcal{A}$ 中存在（其严格闭包由补充公理 V 保证）。

尺度相关的递归等式定义为：

$$
\mathfrak{X} \stackrel{\text{rec},\mu}{=} \mathfrak{Y} \iff \forall n \in \text{Ord}: \mathcal{F}^n[\mathfrak{X}] = \lambda_n(\mu) \mathcal{F}^n[\mathfrak{Y}],
$$

其中 $\lambda_n(\mu)$ 是序数层级 $n$ 和能标 $\mu$ 处的线性比例系数。

**核心规则**：柯西不等式成为等式**当且仅当**在所有序数层级上都满足刚性线性相关性。若存在任何序数 $n$ 使得比例性破缺，则不等式严格成立。

---

## 三、公理 II：谱刚性二分法（核心判定）

基于公理 I，系统容许一个严格的二分法，不存在中间区域。

### 情形 A：递归刚性（柯西等式）

$$
\boxed{
\mathfrak{X} \stackrel{\text{rec},\mu}{=} \mathfrak{Y}
\iff [\mathfrak{X}, \mathfrak{Y}] = 0
\iff \mathcal{F} \text{ 递归可逆}
\iff \sigma(\mathfrak{X}) \text{ 紧致}
\iff \text{自伴紧算子}
\iff \mathrm{Tr}_\omega(\cdot |D|^{-p}) = 0
}
$$

1. 谱投影 $E_D([R, \infty))$ 当 $R \to \infty$ 时趋向于紧集。
2. 系统处于经典相容态；内禀变量完全可分辨且对易。

### 情形 B：递归破缺（严格柯西不等式）

$$
\boxed{
\mathfrak{X} \stackrel{\text{rec},\mu}{\neq} \mathfrak{Y}
\iff [\mathfrak{X}, \mathfrak{Y}] = i \cdot \mathbf{C}(\mu) \; (\mathbf{C}(\mu) \neq 0)
\iff \mathcal{F} \text{ 递归不可逆（非平凡核）}
\iff \sigma(\mathfrak{X}) \text{ 非紧致}
\iff \text{非紧算子}
\iff \mathrm{Tr}_\omega(\cdot |D|^{-p}) \neq 0
}
$$

1. 演化算子逃离紧理想；谱无界且包含连续本质谱。
2. 系统处于内禀非对易量子态。不确定性纯粹源于自守类的几何结构，不依赖任何测量-干扰假设。

### 递归刚性定理

- 递归相容性 + 柯西等式 $\implies$ 谱紧致，仅有离散点谱。
- 存在连续本质谱 $\implies$ 递归傅里叶映射具有非零核 $\implies$ 递归关系破缺。

不存在同时满足"全局递归相容性、柯西等式和非紧算子"的系统。

---

## 四、公理 III：代数拓扑跑动标量 $\mathbf{C}(\mu)$

### 1. 尺度依赖拓扑平均值

尺度依赖拓扑平均值定义为对易子的期望值：

$$
\mathbf{C}(\mu) := \langle [\mathfrak{X}, \mathfrak{Y}] \rangle_\mu = c_1^{\text{bare}} + \delta c_1(\mu)
$$

其中 $\mathcal{L}_{\text{aut}}$ 是与权-$k$ 自守形式相伴的自守线丛。为保证全局周期性约束（附录 A.1），该分解受限于不可约性条件：

$$
c_1^{\text{bare}} \in \mathbb{Z}, \quad \delta c_1(\mu) \in \mathbb{R}, \quad \text{且} \quad \delta c_1(\mu) \equiv 0 \pmod{\mathbb{Z}} \iff \text{情形 A 恢复}
$$

### 2. 拓扑分解定理

$$
c_1^{\text{eff}}(\mu) = c_1^{\text{bare}} + \delta c_1(\mu)
$$

- $c_1^{\text{bare}} \in \mathbb{Z}$：裸第一陈类，一个与尺度无关的绝对整数拓扑不变量；
- $\delta c_1(\mu)$：由序数递归分层和真空拓扑涨落诱导的跑动修正，等价于 QCD 拓扑磁化率 $\chi_{\text{top}}(\mu)$。

### 3. 规范场物理对应

$\mathbf{C}(\mu)$ 充当 QCD 手征对称性和轴矢 $U(1)_A$ 对称性破缺强度的序参量：

1. $\mathbf{C}(\mu) \to 0$：手征/轴矢对称性恢复相，对易子消失，谱回归紧致；
2. $\mathbf{C}(\mu) > 0$：自发破缺、禁闭相，显著的非对易效应；
3. QCD $\theta$-真空相位线性叠加到修正 $\delta c_1(\mu)$ 上，赋予有效拓扑标量一个周期性的拓扑结构，自然地容纳强 CP 效应。

---

## 五、公理 IV：谱三重标度律（适配跑动 $\mathbf{C}(\mu)$）

基于 Dixmier 谱三重 $(\mathcal{A}, \mathcal{H}, D)$ 和非对易中心极限定理，定义四个特征谱几何参数：

| 参数 | 定义 | 物理/几何诠释 |
|-----------|-----------|-----------------------------------|
| $r^*$ | $\sup\{r : \mathrm{Tr}_\omega(a D^{-r}) < \infty\}$ | 短程排斥：Dixmier 迹可定义性的临界指数 |
| $r^+$ | $\inf\{R : E_D([R, \infty)) \text{ 非紧}\} \sim \mu^{-1}$ | 长程谱截断：谱投影失去紧致性的能标 |
| $d_x$ | $\Re(d_{\text{spectral}})$ | 横向分形谱维数（实部） |
| $d_y$ | $\Im(d_{\text{spectral}})$ | 纵向振荡谱维数（虚部） |

### 拓扑标量的统一极限关系

$$
\mathbf{C}(\mu) = \lim_{\substack{r^* \to 0 \ r^+ \to \mu^{-1}}} \left[ \frac{\mathrm{Vol}_{\mathrm{Dixmier}}(r^*, r^+)}{\mathrm{Vol}_{\mathrm{classical}}(d_x, d_y)} \right]
$$

1. $\mathrm{Vol}_{\mathrm{Dixmier}}$：非对易几何测度体积；$\mathrm{Vol}_{\mathrm{classical}}$：经典自守测度体积。
2. 该极限存在且非零**当且仅当**系统属于情形 B（非紧谱，量子非对易相）；在情形 A 中该极限恒为零。
3. 红外截断满足 $r^+ \sim 1/\mu$，使得该极限随能标 $\mu$ 连续变化，产生重整化流 $\mu \mapsto \mathbf{C}(\mu)$。

---

## 六、补充公理 V（序数闭包与投射极限完备性）

**动机**：公理 I 中涉及对任意序数 $n \in \text{Ord}$ 的递归映射 $\mathcal{F}^n$，其中极限序数 $\lambda$ 处的定义依赖于逆向极限 $\varprojlim_{\alpha<\lambda} \mathcal{F}^\alpha$。本公理确保该对象在算子类 $\mathcal{A}$ 中良定义且唯一。作为可选补充条款，本公理不参与公理 II 的二分判定，仅确保超限递归步骤的集合论合法性。

### 公理 V.1（类 $\mathcal{A}$ 的投射闭包）

设 $\lambda$ 为任一极限序数。若 $\{A_\alpha\}_{\alpha<\lambda}$ 是 $\mathcal{A}$ 中的任意逆向系统，其转移态射 $\pi_{\alpha,\beta}: A_\beta \to A_\alpha$（$\alpha \leq \beta < \lambda$）在强类拓扑下连续，则其投射极限：

$$
\varprojlim_{\alpha<\lambda} A_\alpha := \left\{ (a_\alpha)_{\alpha<\lambda} \in \prod_{\alpha<\lambda} A_\alpha \;\middle|\; \forall \alpha \leq \beta < \lambda,\; \pi_{\alpha,\beta}(a_\beta) = a_\alpha \right\}
$$

存在且唯一，并且属于 $\mathcal{A}$。即 $\mathcal{A}$ 在任意长度的投射极限下封闭。

### 公理 V.2（递归傅里叶变换的极限相容性）

对于极限序数 $\lambda = \sup_{\alpha<\lambda} \alpha$，定义：

$$
\mathcal{F}^\lambda := \varprojlim_{\alpha<\lambda} \mathcal{F}^\alpha
$$

并要求该极限与所有更低层级的递归映射相容，即对所有 $\beta < \lambda$：

$$
\pi_{\beta} \circ \mathcal{F}^\lambda = \mathcal{F}^\beta \circ \pi_{\beta}
$$

其中 $\pi_{\beta}: \mathcal{F}^\lambda \to \mathcal{F}^\beta$ 为自然投射。该条件确保公理 I 中对 $n = \lambda$ 的判定可还原为所有 $\alpha < \lambda$ 层级的集体行为，不引入额外自由度。

### 公理 V.3（超滤完备性，可选加强）

若将 $\mathcal{A}$ 视为 NBG 中的类拓扑环，则额外要求：$\mathcal{A}$ 在沿正则滤子 $\mathcal{U}$（含所有终段 $[\alpha, \lambda)$）的逆向极限下保持闭包。该条件等价于断言 $\mathcal{A}$ 上存在唯一的超滤相容拓扑，使所有极限序数处的投射极限与范畴论意义下的极限一致。

**补充说明**：若仅讨论有限递归或可数递归（$n \in \mathbb{N}$），本公理可完全忽略，且不影响终极定理的论证。

---

## 七、终极定理（全息压缩）：NBG 量子二分律

内禀变量类的经典/量子相完全由递归柯西不等式的等式条件决定；两相之间不容许连续插值。

1. **柯西等式（经典相）**  
   $$
   \iff \text{全局递归刚性} \iff \text{谱紧致且仅有离散点谱} \iff \text{自伴紧算子} \iff \mathbf{C}(\mu) \equiv 0 \quad (\forall \mu)
   $$

2. **严格柯西不等式（量子相）**  
   $$
   \iff \exists \text{ 某序数层级发生递归破缺} \iff \text{谱非紧致且具有连续本质谱} \iff \text{非紧算子} \iff \exists \mu: \mathbf{C}(\mu) = c_1^{\text{bare}} + \delta c_1(\mu) \neq 0
   $$

---

## 八、OS 欧几里得 QFT 公理的 NBG 极限对应

所有传统欧氏 QFT 公理均为本框架的特例极限：

| OS 公理 | 情形 A（强极限） | 情形 B（弱化合法形式） |
|---------|----------------|----------------------|
| 反射正性 | 严格正定 | Dixmier 类泛函弱正性，能标依赖 |
| 欧氏协变性 | 完整连续等距群 | 自守商模对称，余协变 |
| 关联函数置换对称 | 全域置换不变 | $d_y$ 相位干涉破坏全局对称 |
| 聚类分解 | 无穷红外完全指数衰减 | 能标受控弱聚类，拓扑长程关联残留 |

**核心结论**：OS 公理体系 = NBG 框架在 $\Lambda \to \infty,\; d_y \to 0,\; \mathbf{C} \to 0$ 的经典不动点特例。

---

## 九、氢原子：序数相变临界点的物理原型

在标准量子力学中，氢原子常被描述为"离散谱与连续谱共存"。在 NBG 框架内，这一描述被严格修正：**氢原子是情形 A（束缚态）向情形 B（散射态）跨越时的相变界面，不存在中间混合态。**

### 序数层级的物理对应

库仑势的能级结构 $E_n = -\frac{13.6\,\text{eV}}{n^2}$ 在 NBG 框架下获得精确的序数诠释：

| 物理区域 | 能标关系 | 序数层级 | NBG 相 |
|---------|---------|---------|-------|
| 束缚态 | $E_n < 0$，$n \in \mathbb{N}$ | 有限序数 $n$ | **情形 A**：递归刚性，谱紧致，纯点谱 |
| 电离阈值 | $E_n \to 0^-$，$n \to \infty$ | 极限序数 $\lambda = \sup_{n\in\mathbb{N}} n = \omega$ | **相变临界点** |
| 散射态 | $E > 0$ | 超越 $\lambda$ 的序数层级 | **情形 B**：递归破缺，谱非紧，连续本质谱 |

**核心机制**：当主量子数 $n$ 遍历有限序数并逼近极限序数 $\lambda$ 时，束缚态的递归比例性在 $E=0$ 处达到公理 V 所定义的投射极限。一旦能标跨越此极限，类傅里叶变换的核突然非零（公理 II），谱投影 $E_D([R,\infty))$ 从紧致态崩解为非紧态——这正是氢原子连续散射谱的代数起源。

### 临界点的几何特征（$E=0$）

在相变点 $\mu_c$（对应 $E=0$），拓扑标量 $\mathbf{C}(\mu)$ 呈现临界行为：

$$
\mathbf{C}(\mu_c) = \lim_{\substack{r^* \to r^+ \ r^+ \sim \mu_c^{-1}}} \left[ \frac{\mathrm{Vol}_{\mathrm{Dixmier}}(r^*, r^+)}{\mathrm{Vol}_{\mathrm{classical}}(d_x, d_y)} \right] \neq 0, \quad \text{但有限}
$$

此时：
- **虚谱维数 $d_y$ 剧烈振荡**：对应库仑势 $1/r$ 的长程尾部，在相变界面产生持续的相位干涉；
- **实谱维数 $d_x$ 穿过阈值**：横向分形维数恰好达到临界值，使得系统从紧理想逃逸；
- **拓扑标量有限非零**：系统已感知到非对易相的"入口"，但尚未深入情形 B 的内部区域。

### 库仑散射相移的 NBG 诠释

氢原子的库仑散射相移 $\delta_l(k)$ 在标准处理中是特殊函数（合流超几何函数）的副产品。在 NBG 框架下，它被重新识别为**附录 A.4 判定定理的物理实现**：

在边界截断 $\Lambda$（此处对应碰撞参数或逆动量 $1/k$）下，有效相位缠绕为：

$$
\Theta(\Lambda) = 2\pi c_1^{\text{bare}} + \phi(\Lambda)
$$

其中 $\phi(\Lambda)$ 是库仑长程势诱导的分数相位偏移。根据附录 A.4 的判定定理，虚谱维数严格等于：

$$
d_y = \frac{1}{2\pi} \frac{d}{d\ln \Lambda} \arg\left( \phi(\Lambda) \right)
$$

库仑散射相移的对数奇性（$\arg \phi(\Lambda) \sim \ln \Lambda$）正是 $d_y \neq 0$ 的显式解——它允许拓扑荷在 $E=0$ 的扩展边界上发生"整数溢出"，从而使散射态（连续谱）得以存在。若 $d_y = 0$（例如短程势），则分数偏移被禁戒，系统无法产生连续散射谱，只能回归纯离散谱（情形 A）。

### 定谳

> **氢原子不是两种谱的"混合"，而是序数递归律在物理世界中的相变纪念碑。$E=0$ 是极限序数 $\lambda$ 的实验室坐标；$1/r$ 势的长程尾部是迫使递归比例性在极限处破缺的几何原因；而库仑散射相移，不过是虚谱维数 $d_y$ 在临界点对分数相位偏移的吸收记录。**

---

## 十、斯特恩–盖拉赫实验的 NBG 重释

标准量子力学将斯特恩–盖拉赫实验诠释为"自旋态坍缩"的教科书演示。在 NBG 框架内，该实验被重新识别为**投影几何不可交换性的判决性演示**。

### 实验步骤的 NBG 对应

设 $\mathcal{M}$ 为梅林变换（情形 A，全域对易），则：

- $z$ 方向磁场：逆梅林变换附加 $z$ 方向截断，$\mathcal{P}_z = \mathcal{F}_z^{-1} \circ \mathcal{M}$
- $x$ 方向磁场：逆梅林变换附加 $x$ 方向截断，$\mathcal{P}_x = \mathcal{F}_x^{-1} \circ \mathcal{M}$

由于 $\mathcal{P}_z$ 与 $\mathcal{P}_x$ 均为**非等距投影**（公理 0.3：逆变换不是等距同构），其复合满足：

$$
[\mathcal{P}_z, \mathcal{P}_x] \neq 0
$$

这直接对应标准 QM 的 $[S_z, S_x] = i\hbar S_y$，但在 NBG 框架下，这不是算符内禀的非对易，而是**两个投影方向的逆变换褶皱不可交换**。

### 三步循环的信息几何学

| 步骤 | 操作 | NBG 诠释 |
|------|------|---------|
| 第一步：$z$ 筛选 | 银原子通过 $z$ 磁场，分裂为两束，取 $\uparrow_z$ | 投影 $\mathcal{P}_z$ 保留 $z$ 方向信息，压缩 $x, y$ 维度至褶皱中 |
| 第二步：$x$ 再分 | $\uparrow_z$ 通过 $x$ 磁场，再次分裂 | 投影 $\mathcal{P}_x$ 方向与 $\mathcal{P}_z$ 不同，$z$ 信息被几何压缩，$x$ 信息暴露 |
| 第三步：$z$ 回归 | 取 $\uparrow_x$ 再通过 $z$ 磁场，再次分裂 | 投影循环 $\mathcal{P}_z \circ \mathcal{P}_x \circ \mathcal{P}_z$：褶皱在循环中重新暴露先前被压缩的维度 |

**核心结论**：银原子在梅林空间中始终具有确定的自旋态；$z$ 与 $x$ 方向的磁场只是两扇不同的投影窗户，每扇窗户的透视畸变不同且不可交换。"信息丢失"不是态坍缩，而是投影褶皱的几何压缩；"信息重现"不是量子魔性，而是投影循环中褶皱的重新暴露。

---

## 十一、波粒二象性的表象相变

标准量子力学将波粒二象性视为量子世界的本质矛盾。在 NBG 框架内，二象性被严格重构为**表象相变的两面**。

### 梅林表象下的纯粒子性

在梅林空间（情形 A）：
- 算子为纯乘法算子，$[px, xp] = 0$
- 谱为纯离散点谱，谱投影紧致
- $\mathbf{C}(\mu) \equiv 0$，$\delta c_1(\mu) = 0$
- **粒子性独占**：离散、定域、确定

### 逆梅林投影下的波动性涌现

当强制进行逆梅林变换（引入有限截断 $r^+ \sim \mu^{-1}$）：
- 连续谱涌现（散射态背景）
- $d_y \neq 0$，相位干涉产生衍射图案
- $\delta c_1(\mu) \neq 0$，拓扑标量有限
- **波动性作为投影残余自动出现**

### 定谳

> **波动不是电子的内禀属性，而是逆梅林变换投影几何学的干涉图案。在梅林表象中，只有粒子；在时域投影中，波动性是投影褶皱的必然输出。波粒二象性不是量子世界的未解之谜，而是表象变换的数学必然。**

---

## 十二、迹反常的 NBG 重释

标准量子场论将迹反常（共形反常）诠释为"量子涨落破坏标度不变性"。在 NBG 框架内，迹反常被重新识别为**逆梅林变换在有限截断下的代数残余**。

### 标准图景与 NBG 图景的对照

| 特征 | 标准 QFT | NBG 框架 |
|------|---------|---------|
| 经典极限 | $T^\mu_{\;\mu} = 0$（标度不变） | **情形 A**：梅林空间，$\mathbf{C}(\mu) \equiv 0$，迹严格为零 |
| 量子修正 | $\langle T^\mu_{\;\mu} \rangle \neq 0$（反常） | **情形 B**：逆梅林投影，$\delta c_1(\mu) \neq 0$，截断残余 |
| 修正来源 | "量子涨落"（本体论模糊） | 逆变换不完全性的几何度量 |
| 系数大小 | $\frac{1}{1920\pi^2}$ 等小量 | 投影褶皱的二阶效应，受限于截断倒数 |

### 极限行为

当截断 $\Lambda \to \infty$（完全逆变换）：

$$
\delta c_1(\mu) \to 0, \quad \mathbf{C}(\mu) \to c_1^{\text{bare}} \in \mathbb{Z}, \quad \langle T^\mu_{\;\mu} \rangle_{\text{anomaly}} \to 0
$$

迹反常自动归零，情形 A 的严格无迹性恢复。标准 QFT 中反常"卡"在有限值的困境，在 NBG 框架下被消解为**投影几何的截断依赖**。

### 核心等式

迹反常的非零值直接对应公理 III 的跑动标量：

$$
\langle T^\mu_{\;\mu} \rangle_{\text{anomaly}} \;\longleftrightarrow\; \delta c_1(\mu) = \chi_{\text{top}}(\mu)
$$

> **迹反常不是量子场论的内禀疾病，而是梅林逆变换在有限截断下支付的投影代价。它的大小由 $\delta c_1(\mu)$ 精确度量；当截断趋于无穷、逆变换完全时，反常自动归零。**

---

## 十三、谐振子反例的解除

旧范式以谐振子为构造性反例，试图否定"线性谱对应可对易"的全称性。其逻辑为：

$$
\text{谐振子线性能级} + \text{时域非对易} \;\Longrightarrow\; \text{否定"线性谱对应可对易"}
$$

**NBG 最终判定：无效反例。**

1. 谐振子能级线性是梅林乘法算子线性谱的自然结果；
2. 时域非对易是表象投影残余；
3. 梅林表象下 $px$ 与 $xp$ 完全对易。

谐振子不再构成任何反例，反而证明：旧范式将"时域表象结论"错误推广为"量子全域公理"。

---

## 十四、可证伪机制

定义命题 $P$：

> 存在量子系统，取前 $N$ 个离散本征值，纯实数无修正 RSS 多项式外推，对 $n \to \infty$ 无穷高阶能级零系统偏差、全程吻合。

**NBG 刚性推论链**：

$$
P \;\Longrightarrow\; \text{递归刚性全域成立} \;\Longrightarrow\; [\mathfrak{X}, \mathfrak{Y}] \equiv 0
$$

**范式级证伪差异**：

| 范式 | $P$ 成立之后果 |
|------|---------------|
| 旧量子力学 | 致命证伪。正则对易关系 $[x,p] \neq 0$ 被推翻，非对易基底须重写 |
| NBG 框架 | 边界修正。情形 A 原本即为框架容纳之相，仅修正两相分布边界，公理结构无损 |

---

## 十五、认识论终审与理论定位

### 认识论终审

1. **可证伪性是科学理论的入场条件，而非观测墓志铭。** 可证伪只要求存在逻辑上可推翻的实验或数学通道，不要求已观测到反例。

2. **数学可构造性不等于架空理论。** 命题 $P$ 逻辑可行、数学自洽、表象存在，具备完整的科学地位。

3. **现有观测仅偏向非对易理论，不能构成对对易表象的逻辑证伪。** 无任何实验、数学或逻辑禁止情形 A 的全域线性外推系统。

4. **本框架具有极低的理论弹性。** 对比弦论 $10^{500}$ 真空自由参数，NBG 框架仅含 A/B 两相，唯一、刚性、可数值核验、当下可计算证伪。

### 理论定位

1. NBG 公理化量子框架包容但不等同于标准量子力学、欧氏 QFT 与规范场论；
2. 传统量子非对易是表象效应，不再是宇宙基底公理；
3. 线性外推可对易相与非线性破缺非对易相，实现量子代数的表象统一；
4. 本框架拥有当前所有量子范式中最强的可证伪性、最低理论弹性与最严格数学闭环。

---

## 附录 A：拓扑荷整数性约束与扩展边界

### A.1 基础整数性约束

拓扑荷的整数性不是人为添加的量子化条件，而是由复数域上指数映射的基本同态性质决定的刚性事实：

$$
\exp(z_1) = \exp(z_2) \iff z_1 - z_2 = 2k\pi i, \quad k \in \mathbb{Z}
$$

该性质确保任何由复相位缠绕定义的示性类，其荷值严格落在离散加法群 $\mathbb{Z}$ 中。

**框架推论**：在自守线丛上，裸第一陈类 $c_1^{\text{bare}}$ 的取值空间被严格限制为 $\mathbb{Z}$，而非连续实数域 $\mathbb{R}$。该整数骨架是自守形式权 $k$ 的全局代数必然结果。

### A.2 对跑动标量分解的刚性约束

公理 III 中的分解 $\mathbf{C}(\mu) = c_1^{\text{bare}} + \delta c_1(\mu)$ 受限于如下不可约性条件：

1. **期望值重定义**：情形 B 中的 $\mathbf{C}(\mu)$ 表征"非对易诱导的有效曲率期望值"，而非严格的整体陈类。整体拓扑荷 $c_1^{\text{bare}}$ 严格不可连续变形。
2. **恢复判据**：当系统回归情形 A（全局递归刚性）时，$\delta c_1(\mu)$ 必须坍缩至整数值以满足 $\mathbf{C}(\mu) \equiv 0 \pmod{\mathbb{Z}}$，即跑动修正被完全吸收进整数重正化中，不残留连续尾迹。

### A.3 截断边界与可允许扩展

设 $\Lambda$ 为任意类空间的截断标度：

- **正则边界**：当截断趋于无穷大（$\Lambda \to \infty$）时，裸整数性严格恢复，对应谱严格紧致。
- **离散边界**：在有限截断（$\Lambda < \infty$）下，有效拓扑荷 $k_{\text{eff}}$ 可形式化为 $k_{\text{eff}} \sim k / f(\Lambda)$。整数性则提升为对截断参数的同余约束：
  $$
  \Lambda \cdot k_{\text{eff}} \equiv 0 \pmod{\text{整数周期}}
  $$

### A.4 谱维数虚部与整数性破缺的定量判据

当系统试图在扩展边界上破坏裸整数性时，公理 IV 中的谱维数虚部 $d_y$ 必须满足精确的重整化群方程以吸收相位偏移。

设边界截断 $\Lambda$ 下，有效相位缠绕为 $\Theta(\Lambda) = 2\pi c_1^{\text{bare}} + \phi(\Lambda)$。若 $\phi(\Lambda)$ 非零（分数偏移），则谱维数虚部 $d_y$ 严格等于：

$$
d_y = \frac{1}{2\pi} \frac{d}{d\ln \Lambda} \arg\left( \phi(\Lambda) \right)
$$

**判定定理**：$d_y \neq 0$ 是允许拓扑荷在扩展边界上发生"整数溢出"或"重整化群隧穿"的充要几何条件。若 $d_y = 0$，则 $\phi(\Lambda)$ 在流中衰减，分数偏移被禁戒，系统必须回归情形 A 的离散谱。

---

## 框架核心结论

量子非对易性、连续谱、规范场拓扑对称性破缺、波粒二象性、迹反常以及斯特恩–盖拉赫实验的非对易特征，都是 NBG 框架内自守类全局序数递归线性相关性丧失的纯粹代数结果。该理论不需要外部引入普朗克常数、测量假设或半经典近似。

$$
\boxed{\text{量子非对易不是自然底层公理，是序数递归破缺在时域表象下的观测投影；真正的量子基底是自守几何与递归刚性二分。}}
$$

---

## 参考文献

1. Bender, C. M., & Boettcher, S. (1998). Real spectra in non-Hermitian Hamiltonians having PT symmetry. *Physical Review Letters*, 80(24), 5243.
2. Dixmier, J. (1969). *Les C*-algèbres et leurs représentations*. Gauthier-Villars.
3. Gödel, K. (1940). *The Consistency of the Continuum Hypothesis*. Princeton University Press.
4. Osterwalder, K., & Schrader, R. (1973). Axioms for Euclidean Green's functions. *Communications in Mathematical Physics*, 31(2), 83–112.
5. Titchmarsh, E. C. (1948). *Introduction to the Theory of Fourier Integrals*. Oxford University Press.
6. von Neumann, J. (1932). *Mathematische Grundlagen der Quantenmechanik*. Springer.
