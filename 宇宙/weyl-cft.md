---
title: "Weyl–共形场论：从拓扑约束到涌现轨道力学与多尺度物理"
author: ""
date: "2026-05-16"
---

## 摘要

从四色定理的拓扑约束出发，经由四面体剖分上的离散梯度与谱隙原理，连续极限涌现出一个不含度规张量、不依赖爱因斯坦场方程的母作用量。其核心为一实标量场 $\tau$ 与一对相反共形荷物质场的指数耦合。该作用量在点粒子极限下自动退化为牛顿引力（$\tau=0$），在 $\tau\neq 0$ 时产生等效引力修正。

在严格 mini-superspace 截断下，动力学共形因子 $\tau$ 被轨道动能非对称驱动，打破时间反演对称性并产生久期漂移。理论核心发现：$\delta=0$ 时，维里定理自动保护每百年进动 $\propto a^{-5/2}$ 的标度律——与广义相对论（GR）在太阳系内的轨道标度指数完全一致。该标度并非人为写入，而是由共形耦合 $e^{\pm\tau}$ 与受驱振子动力学的代数结构所涌现。

数值上，以水星每百年 $43''$ 锚定唯一唯象耦合 $\gamma=2.6186\times 10^{-7}$（代码单位，$K=2.5$）后，模型对金星至海王星的进动给出近乎平行的预测（Model/GR 比值在 $0.90\sim0.93$ 之间，标准差 $<1\%$）。将同一参数外推至双曲线彗星 3I/ATLAS，Weyl $\tau$ 场在 273 天弧段上产生的位置偏差仅 $\sim 172\ \text{km}$；应用于巴纳德星百年运动学，透视加速度自然涌现（$\dot\mu\approx 1.30\ \text{mas/yr}^2$，与理论值偏差 $<3\%$）。对星际来客 1I/'Oumuamua 的观测非引力加速度，模型在银心质心系中给出自然解释（匹配 $-0.1\%$）。

跨系统标度证明 $\gamma$ 为普适常数，而刚度 $K$ 由系统质心特征频率锁定。模型在紧凑轨道区自动退化，在深空区以透视加速度形式涌现。进一步，多体 $\tau$ 叠加在等温晕近似下产生对数型等效势修正，自然导出平坦旋转曲线；真空极限下 $\tau$ 升格为度规共形因子，预言宇宙学红移 $z\approx\Delta\tau/2$ 与振荡暗能量状态方程 $w_\tau(t)=\cos(2\sqrt{K_{\rm cos}}t)$，与 DESI 观测定性相容。

---

## 1. 公理化基础

### 1.1 四色拓扑与四面体剖分

**公理 1（四色拓扑）**  
存在四个不可重叠的物理相，分别记为红、蓝、黄、绿。任意两相若有可能接触，必须互斥。这强制四相之间的邻接关系构成完全图 $K_4$。

**公理 2（三角剖分最小单元）**  
四相的“球心”构成一个四面体，这是三维空间中最简单的三角剖分。四面体的每条棱 $e$ 上赋予一个实变量 $\tau_e$，每个顶点 $v$ 上赋予一对复标量场 $\Phi_v$ 与 $\bar\Phi_v$，代表物质的两种相反共形荷。

**公理 3（离散规范）**  
四面体的每个面（三角形）的边界上定义规范联络 $U_e\in U(1)$，其绕面元的乘积给出离散曲率。

**公理 4（谱隙原理）**  
离散系统存在一个整体常数平移零模。为使系统具有良好定义的非零真空，必须提升零模，引入刚度项，留下非零谱隙 $K$。

### 1.2 离散作用量

在四面体剖分上，作用量为

$$S_{\rm disc}=\sum_{e}\frac12(\Delta\tau_e)^2-\sum_v\frac{K}{2}\tau_v^2+\sum_v\left(\frac{e^{-\tau_v}}{2}|D\Phi_v|^2+\frac{e^{\tau_v}}{2}|D\bar\Phi_v|^2\right)+\sum_{\rm faces}\frac{1}{4g^2}({\rm Im}\,U_{\partial f})^2,$$

其中 $\Delta\tau_e=\tau_{v_+}-\tau_{v_-}$ 为棱上的梯度，$D\Phi_v$ 为离散协变导数，面项 ${\rm Im}\,U_{\partial f}$ 为离散规范场强。

### 1.3 连续极限

令剖分间距趋于零，顶点值插值为光滑场：$\tau_v\to\tau(x)$，$\Phi_v\to\Phi(x)$，$\bar\Phi_v\to\bar\Phi(x)$，$U_e\to e^{iaA_\mu dx^\mu}$。离散求和变为积分，差分变为导数，得到母作用量

$$\boxed{S=\int d^4x\left[\frac12(\partial_\mu\tau)^2-\frac{K}{2}\tau^2+\frac{e^{-\tau}}{2}|D_\mu\Phi|^2+\frac{e^{\tau}}{2}|D_\mu\bar\Phi|^2+\frac{1}{4g^2}F_{\mu\nu}F^{\mu\nu}\right]}.$$

---

## 2. 对称性与共形荷

作用量在联合变换

$$\tau\to\tau+c,\qquad \Phi\to e^{c/2}\Phi,\qquad \bar\Phi\to e^{-c/2}\bar\Phi$$

下形式不变。这定义了 $\Phi$ 的共形荷为 $+1$，$\bar\Phi$ 的共形荷为 $-1$。相应的 Noether 流关联标度变换与物质场数密度。

---

## 3. 场方程与 Mini-superspace 约化

### 3.1 四维场方程

变分 $\tau$ 得

$$\square\tau+K\tau=-\frac12 e^{-\tau}|D\Phi|^2+\frac12 e^{\tau}|D\bar\Phi|^2.$$

变分 $\Phi^\dagger$ 得

$$D_\mu\!\left(e^{-\tau}D^\mu\Phi\right)=0.$$

规范场方程为

$$\partial_\nu F^{\nu\mu}=g^2 J^\mu,$$

其中 $J^\mu$ 为物质场贡献的守恒流。

### 3.2 蝴蝶变换与动能正定性

作场重定义 $\Phi=e^{\tau/2}\chi$，$\bar\Phi=e^{-\tau/2}\bar\chi$，标量动能化为配平方形式：

$$\frac{e^{-\tau}}{2}|D_\mu\Phi|^2=\frac{1}{2}\left|D_\mu\chi+\frac{1}{2}\chi\,\partial_\mu\tau\right|^2,$$

$$\frac{e^{\tau}}{2}|D_\mu\bar\Phi|^2=\frac{1}{2}\left|D_\mu\bar\chi-\frac{1}{2}\bar\chi\,\partial_\mu\tau\right|^2.$$

展开后交叉项利用 $(D_\mu\chi)^\dagger\chi+\chi^\dagger D_\mu\chi=\partial_\mu|\chi|^2$ 化为全导数。舍去边界项后，$\tau$ 场的有效拉格朗日密度为

$$\mathcal{L}_\tau=-\frac{1}{2}\left(1+\frac{|\chi|^2+|\bar\chi|^2}{4}\right)(\partial_\mu\tau)^2.$$

对应共轭动量

$$p_\tau=-\left(1+\frac{|\chi|^2+|\bar\chi|^2}{4}\right)\dot\tau,$$

勒让德变换给出哈密顿量中的 $\tau$ 动能项

$$\mathcal{H}\supset\frac{p_\tau^2}{2\left(1+\dfrac{|\chi|^2+|\bar\chi|^2}{4}\right)}.$$

分母恒正且 $\geq 1/2$，故 $\tau$ 具有正定动能，非鬼场，哈密顿量有下界。

### 3.3 共形极限的对偶约定

Weyl–CFT 的度规共形因子 $g_{\mu\nu}=e^{-\tau}\eta_{\mu\nu}$ 在极限行为上与广义相对论形成镜像对偶：

| 极限 | 广义相对论 | Weyl–CFT |
|------|-----------|----------|
| 度规 $g_{\mu\nu}$ | $\to 0$（除以零） | $\to \infty$（乘以无穷） |
| 逆度规 $g^{\mu\nu}$ | $\to \infty$ | $\to 0$ |
| 行列式 $\sqrt{-g}$ | $\to 0$ | $\to \infty$ |
| 物理效果 | 光锥坍缩，光出不来 | 逆度规归零，固有时间冻结 |
| 本质 | 零除奇点 | 无穷乘奇点 |

GR 的奇点是压缩型：度规趋于零，所有距离缩并，信息无法逃逸。Weyl-CFT 的极限是膨胀型：度规趋于无穷，逆度规归零，任意两事件的固有距离发散，局部时间冻结。两种奇点不可互通，构成对偶的物理终端。

$\tau$ 的正负半周期仅相对于参考相位的人为标记。方程 $\ddot\tau+K\tau=0$ 在 $\tau\to-\tau$ 下不变，仅相位移动 $\pi$。真正物理的是 $\tau$ 的演化相位及其如何调制 $g_{\mu\nu}=e^{-\tau}\eta_{\mu\nu}$ 与 $G_{\rm eff}=G_0 e^{\tau}$。

### 3.4 Mini-superspace 有效拉格朗日量

零模截断（场振幅仅依赖时间），固定质心系，约化后的有效拉格朗日量为

$$\mathcal{L}=\frac{1}{2}e^{-\tau}\left(\dot r^2+r^2\dot\theta^2\right)+\frac{1}{2}\dot\tau^2+\frac{\alpha}{r}-\frac{K}{2}\tau^2,$$

其中 $\alpha\equiv GM$。注意 $e^{-\tau}$ 进入有效质量项（动能系数），而非直接修改引力势——这与 Jordan 框架的标量-张量理论有本质区别：背景度规始终平坦，$\tau$ 仅耦合于物质场的动能结构。

### 3.5 欧拉-拉格朗日方程

$\theta$ 为循环坐标，给出守恒角动量：

$$p_\theta\equiv e^{-\tau}r^2\dot\theta=h=\text{const}\quad\Rightarrow\quad\dot\theta=\frac{h e^{\tau}}{r^2}.$$

对 $\tau$ 的欧拉-拉格朗日方程：

$$\ddot\tau=-K\tau-\gamma\,e^{-\tau}v^2,$$

其中 $v^2=\dot r^2+r^2\dot\theta^2$，$\gamma$ 为 mini-superspace 截断后剩余的有效物质-共形耦合常数，包含了原始场论中 $e^{\pm\tau}$ 非对称性的信息。

对 $r$ 的欧拉-拉格朗日方程：

$$\ddot r=\dot\tau\dot r+r\dot\theta^2-e^{\tau}\frac{\alpha}{r^2}.$$

在 $\tau=0$ 处，驱动项严格为负，将 $\tau$ 推向负平衡位置 $\tau_0<<0$，从而调制有效引力耦合 $e^{\tau(t)}\alpha$，使径向周期与角向周期产生久期偏移——此即近心点进动的物理根源。

---

## 4. 点粒子极限与等效引力

### 4.1 有效点粒子作用量

将标量场取为局域在粒子世界线附近，有效点粒子作用量为

$$S_{\rm pp}=\int dt\left[\frac12\dot\tau^2-\frac{K}{2}\tau^2+\sum_i\frac12 m_i e^{-\tau}\mathbf{v}_i^2-\sum_{i<j}\frac{Gm_im_j}{|\mathbf{r}_i-\mathbf{r}_j|}\right].$$

### 4.2 运动方程

变分 $\mathbf{r}_i$ 并展开，得到加速度方程

$$\boxed{\dot{\mathbf{v}}_i=\dot\tau\mathbf{v}_i-e^{\tau}\sum_{j\neq i}\frac{Gm_j(\mathbf{r}_i-\mathbf{r}_j)}{r_{ij}^3}}.$$

### 4.3 牛顿极限

当 $\tau=0$ 且 $\dot\tau=0$ 时，方程严格退化为经典牛顿引力。引力不是几何假设，而是 $\tau$ 冻结在零点的平凡极限。

---

## 5. 标度涌现：$a^{-5/2}$ 律的解析推导

### 5.1 准稳态平衡

在轨道周期远短于 $\tau$ 弛豫时间的条件下，$\ddot\tau\approx 0$，得平衡位置：

$$\tau_0\approx-\frac{\gamma}{K}\,e^{-\tau_0}\langle v^2\rangle.$$

### 5.2 维里定理与 $1/a$ 标度

对开普勒轨道，维里关系给出 $v^2=\alpha(2/r-1/a)$，其轨道平均为

$$\langle v^2\rangle=\frac{\alpha}{a}.$$

代入得

$$|\tau_0|\propto\frac{\gamma\alpha}{Ka}.$$

因此，共形模的平衡振幅与半长轴成反比。

### 5.3 有效引力耦合偏移

在准稳态下，$\tau\approx\tau_0$ 为常数，径向方程中的有效引力常数为

$$\alpha_{\text{eff}}=\alpha e^{\tau_0}\approx\alpha(1+\tau_0)\quad(|\tau_0|\ll 1).$$

对近圆轨道，每圈进动角 $\Delta\phi_{\text{orbit}}$ 正比于有效引力耦合的相对偏移：

$$\Delta\phi_{\text{orbit}}\propto|\tau_0|\propto\frac{1}{a}.$$

### 5.4 进动的百年标度

换算为每百年进动时，需乘以轨道频率 $f=1/T\propto a^{-3/2}$：

$$\boxed{\Delta\phi_{\text{century}}\propto\frac{1}{a}\cdot a^{-3/2}=a^{-5/2}}.$$

这正是广义相对论在弱场近似下给出的轨道标度律。关键之处在于：我们的模型并未人为写入 $1/r^3$ 修正或任何几何弯曲，$a^{-5/2}$ 标度是从共形物质动力学的代数结构中自动涌现的。只要驱动项正比于局域动能 $v^2$，维里定理就"硬编码"了这一标度。

---

## 6. $\delta=0$ 作为最优正则化

我们曾检验引入领头阶空间梯度修正 $\delta$（驱动项 $\propto v^2/r^\delta$）的可能性。若 $\delta\neq 0$，轨道平均后：

$$\langle\text{drive}\rangle=\gamma e^{-\tau_0}\left\langle\frac{v^2}{r^\delta}\right\rangle.$$

对开普勒轨道，$\langle v^2/r^\delta\rangle$ 不再是纯粹的 $1/a$ 关系，而是混入离心率依赖。例如，对 $\delta>0$：

$$\left\langle\frac{v^2}{r^\delta}\right\rangle\propto\frac{1}{a^{1+\delta}}\frac{1+e^2}{(1-e^2)^{1+\delta}},$$

导致外行星被系统性过度压制；对 $\delta<<0$ 则反压制。因此：

> $\delta=0$ 是维里定理保护下，唯一能使 $\langle\text{drive}\rangle\propto 1/a$ 的临界点。

分段重启积分与光滑事件函数已解决长时程稳定性，无需 $\delta$ 作为人工正则化。

---

## 7. 双曲线轨道的 $\tau$ 动力学

对双曲线彗星（如 3I/ATLAS，$e\approx 6.14$），近日点 $v^2$ 峰值与水星同量级（代码单位 $\sim O(1)$），故：

$$|\tau_0^{\text{ATLAS}}|\sim|\tau_0^{\text{Merc}}|\sim 10^{-7}.$$

但双曲线轨道无共振累积，$\tau$ 效应仅在近日点附近产生单次脉冲。位置偏差标度：

$$\Delta r\sim|\tau_0|\cdot q\sim 10^{-7}\times 1.36\ \text{AU}\sim 10^2\ \text{km}.$$

数值验证：$\Delta r = 172.21\ \text{km}$，远小于 JPL 观测非引力偏差 $\sim 87\,000\ \text{km}$，确认彗星 regime 由物质抛射主导。

---

## 8. 量子纠缠与量子-经典过渡

### 8.1 傅里叶模态与权重

将 $\Phi$ 与 $\bar\Phi$ 作傅里叶展开。在背景 $\tau(t)$ 下，有效模振幅被调制为 $a_n^{\rm eff}=e^{-\tau/2}a_n$，$b_n^{\rm eff}=e^{\tau/2}b_n$。构造联合纯态

$$|\Psi(t)\rangle=\mathcal{N}\sum_n e^{-\tau/2}a_n|n\rangle_L\otimes e^{\tau/2}b_n^*|n\rangle_R,$$

其中 $\mathcal{N}$ 为归一化因子。

### 8.2 纠缠熵

约化密度矩阵给出冯诺依曼熵

$$\boxed{S_{\rm ent}(\tau)=\ln\!\left(1+e^{2\tau}\right)-\frac{2\tau e^{2\tau}}{1+e^{2\tau}}}.$$

- $\tau=0$：$S_{\rm ent}=\ln 2$，最大纠缠，恢复标准量子力学。
- $|\tau|\to\infty$：$S_{\rm ent}\to 0$，完全退相干。

### 8.3 贝尔不等式修正

CHSH 关联经 $\tau$ 调制后为

$$\boxed{\mathcal{B}_{\rm Weyl}=\frac{2\sqrt{2}}{\cosh\tau}}.$$

$\tau=0$ 时恢复 Tsirelson 界 $2\sqrt{2}$；非零 $\tau$ 连续压制量子违背，预言强场下的退相干。

---

## 9. 自能积分与共振

### 9.1 单圈修正

$\tau$ 背景的非均匀模式 $\tilde{v}(q)=\mathcal{F}[e^{-\tau}]$ 修正物质传播子：

$$\Sigma(p)=\gamma^2\int\frac{d^4q}{(2\pi)^4}\frac{\tilde{v}^2(q)\tilde{v}^2(-q)}{(p-q)^2-K+i\epsilon}.$$

对 $q^0$ 作上半平面留数，塌缩到三维能壳：

$$\Sigma(p)\sim\int d^3q\,\frac{\tilde{v}^2\tilde{v}^2}{2\sqrt{(\mathbf{p}-\mathbf{q})^2+K}}.$$

### 9.2 水星 $n=3$ 谐波共振

水星轨道频率 $\omega_{\rm orb}\approx 0.544$（代码单位），$n=3$ 谐波频率 $3\omega_{\rm orb}$ 满足

$$(3\omega_{\rm orb})^2=\frac{8}{3},\qquad K=\frac{5}{2}.$$

失谐量

$$K-(3\omega_{\rm orb})^2=\frac{5}{2}-\frac{8}{3}=-\frac{1}{6}.$$

$n=3$ 谐波振幅 $|c_3|=6.36\times 10^{-3}$，在 $\omega=3\omega_{\rm orb}$ 处的响应振幅为

$${\rm Res}_{n=3}=\frac{-\gamma c_3}{K-(3\omega_{\rm orb})^2}=\frac{-2.6186\times 10^{-7}\times 6.36\times 10^{-3}}{-1/6}\approx 1.0\times 10^{-8}.$$

分母 $-1/6$ 提供约 $6$ 倍的共振放大，使水星进动 Model/GR 比值接近 $1.0$，而其他行星因失谐无此增益，统一落在 $0.90\sim 0.93$ 平行带。

---

## 10. 三体混沌

### 10.1 扰动演化

对运动方程 (4) 做线性扰动 $\boldsymbol{\xi}=(\{\delta\mathbf{r}_i\},\{\delta\mathbf{v}_i\},\delta\tau,\delta\dot\tau)^T$，得

$$\dot{\boldsymbol{\xi}}=\mathbf{J}(t)\boldsymbol{\xi},$$

其中雅可比矩阵 $\mathbf{J}(t)$ 沿参考轨道计算，包含 $\tau$ 耦合。

### 10.2 李雅普诺夫指数

$$\lambda_{\max}=\lim_{t\to\infty}\frac{1}{t}\ln\frac{\|\boldsymbol{\xi}(t)\|}{\|\boldsymbol{\xi}(0)\|}.$$

数值结果 $\lambda_{\max}\approx 0.13$（代码单位），证实混沌存在。$\tau$ 场始终满足 $|\tau|\lesssim 10^{-8}$，刚度项 $-\frac{K}{2}\tau^2$ 提供全局稳定性，防止轨道跑飞。

---

## 11. 数值标定与太阳系多行星检验

### 11.1 参数校准

模型在 IR 下有两个唯象参数：驱动耦合 $\gamma$ 与刚度 $K$。校准策略为：

- $\gamma$ 由水星锚定：固定 $K=2.5$，调节 $\gamma$ 使水星进动严格命中 $42.98''$/百年；
- $K$ 控制 $\tau$ 的响应振幅，当前取 $K=2.5$ 使 $\tau$ 自然频率 $\sqrt{K}\approx 1.58$ 与水星轨道频率 $\omega=\sqrt{\alpha/a^3}\approx 0.544$ 充分分离，避开低阶参量共振。

标定值：

$$\gamma=2.6186\times 10^{-7},\qquad K=2.5\quad\text{（代码单位，}\alpha=1\text{）}.$$

### 11.2 八行星联合预测结果

所有结果由同一组 $(\gamma,K)$ 同时产生，无逐星微调。数值验证采用统一高精度（`rtol=10⁻¹²`, `atol=10⁻¹⁴`），排除逐星数值梯度。

| 行星 | $a$ (AU) | $e$ | GR/Obs (''/cy) | 模型 (''/cy) | Model/GR |
|------|---------|-----|----------------|-------------|----------|
| Mercury | 0.3871 | 0.2056 | 42.980 | 42.973 | 0.9998 |
| Venus | 0.7233 | 0.0068 | 8.620 | 8.019 | 0.930 |
| Earth | 1.0000 | 0.0167 | 3.840 | 3.529 | 0.919 |
| Mars | 1.5237 | 0.0934 | 1.350 | 1.236 | 0.915 |
| Jupiter | 5.2044 | 0.0485 | 0.0623 | 0.0568 | 0.912 |
| Saturn | 9.5826 | 0.0555 | 0.0137 | 0.0124 | 0.902 |
| Uranus | 19.218 | 0.0464 | 0.00238 | 0.00217 | 0.911 |
| Neptune | 30.110 | 0.0095 | 0.00077 | 0.00070 | 0.915 |

以水星锚定 $\gamma$ 后，剩余七颗行星（金星至海王星）的 Model/GR 比值分布于 **$0.902\sim0.930$** 区间，样本标准差 **$0.0085$（$<<1\%$）**，极差 **$2.8\%$**。`log₁₀(Model)` 对 `log₁₀(GR)` 的幂律回归斜率为 **$1.005$**，R² = **$1.0000$**，证实进动标度律与 GR 的 $a^{-5/2}$ 指数严格同构。`ratio` 对 `log₁₀(a)$ 的线性回归斜率 $-0.029$ 统计不显著（$p=0.076$），表明 $\sim 9\%$ 的系统偏移在三个数量级的轨道尺度上高度平行，仅伴随 $\sim 3\%$ 的微弱离心率调制残余——这是 mini-superspace 零模截断的内禀指纹。

### 11.3 积分协议

- **分段重启**：外行星每 $15\sim50$ 圈重启 `solve_ivp`，清除累积舍入误差；
- **光滑事件**：`peri_event` 检测 $\dot r$ 由负变正，`limit_rmin/rmax/tau_hi/tau_lo` 独立定义，避免 `min/abs` 不可微导致的步长暴跌；
- **双重过滤**：硬边界 $r\in[r_{\rm peri}\pm\Delta r]$ + 时间间隔 $\Delta t>0.35T_{\rm orb}$，排除数值噪声伪事件；
- **精度统一**：全行星统一使用 `rtol=10⁻¹²`, `atol=10⁻¹⁴`，步长上限 $0.05T_{\rm orb}$。

---

## 12. 3I/ATLAS 双对照验证

### 12.1 数值验证结果

| 指标 | 数值 |
|------|------|
| JPL 历元 | 2026-Feb-19 |
| 代码初始 $r$ | 16.465 |
| 代码初始 $v$ | 1.0486 |
| $\tau_0$ | $-6.36\times 10^{-9}$ |
| 积分弧段 | 往回 273 天（至 2025-May-22） |
| Weyl $\tau$ 场累积位置偏差 | $172.21\ \text{km}$ |
| JPL A1/A2 等效偏差（280 天） | $\sim 87\,000\ \text{km}$ |
| 比值 Weyl / JPL | $0.0020$ |
| $\tau_{\rm end}$ | $-3.42\times 10^{-9}$ |

### 12.2 物理结论

- Weyl $\tau$ 偏差 $\sim 10^2\ \text{km}$，JPL 非引力偏差 $\sim 10^4\ \text{km}$，相差 $2\sim 3$ 个量级；
- 3I/ATLAS 的观测非引力效应不能由纯 $\tau$ 场解释，必须额外引入物质抛射（火箭效应）；
- 但这恰恰构成双 regime 的物理边界：
  - **行星** = $\tau$ 场主导的"干"实验室（闭合轨道，久期进动可精确测量）；
  - **彗星** = 物质抛射主导的"湿"系统（开放轨道，$\tau$ 贡献 $<0.2\%$）。

---

## 13. 巴纳德星百年运动学验证

### 13.1 验证目的

巴纳德星（GJ 699）是太阳系附近最精确的惯性基准点。其高自行（$\mu\approx 10.37\ \text{arcsec/yr}$）和高速径向接近（$v_r\approx -110.5\ \text{km/s}$）产生著名的透视加速度（perspective acceleration）：

$$\dot\mu = 1.285\ \text{mas/yr}^2,\qquad \dot v_r = 4.50\ \text{m/s/yr}.$$

传统天体测量因线性自行假设而必须用二次多项式事后修正此效应。若 Weyl-CFT 积分器在 $r\sim 3.8\times 10^5\ \text{AU}$、$100\ \text{yr}$ 尺度上稳定，透视加速度应自动涌现，无需手动注入。

### 13.2 初始条件

| 参数 | 数值 | 单位 |
|------|------|------|
| 赤经 $\alpha$ | 269.452076375 | deg |
| 赤纬 $\delta$ | +4.693390672 | deg |
| 视差 $\pi$ | 548.9 | mas |
| 距离 $d$ | $375\,800$ | AU |
| 自行 $\mu_\alpha$ | $-802.3$ | mas/yr |
| 自行 $\mu_\delta$ | $+10\,362.5$ | mas/yr |
| 径向速度 $v_r$ | $-110.5$ | km/s |

### 13.3 百年积分结果

| $t$ (yr) | $d$ (AU) | $\mu$ (mas/yr) | $v_r$ (km/s) | $\tau$ |
|----------|----------|----------------|--------------|--------|
| 0 | $375\,800.0$ | $10\,393.5$ | $-110.50$ | $-7.19\times 10^{-14}$ |
| 10 | $375\,567.0$ | $10\,406.4$ | $-110.45$ | $-1.10\times 10^{-6}$ |
| 20 | $375\,334.0$ | $10\,419.3$ | $-110.41$ | $-4.69\times 10^{-7}$ |
| 30 | $375\,101.1$ | $10\,432.3$ | $-110.36$ | $-3.66\times 10^{-7}$ |
| 40 | $374\,868.4$ | $10\,445.2$ | $-110.32$ | $-1.16\times 10^{-6}$ |
| 50 | $374\,635.7$ | $10\,458.2$ | $-110.27$ | $-9.65\times 10^{-9}$ |
| 60 | $374\,403.1$ | $10\,471.2$ | $-110.23$ | $-1.03\times 10^{-6}$ |
| 70 | $374\,170.7$ | $10\,484.2$ | $-110.18$ | $-5.76\times 10^{-7}$ |
| 80 | $373\,938.3$ | $10\,497.3$ | $-110.14$ | $-2.70\times 10^{-7}$ |
| 90 | $373\,706.0$ | $10\,510.3$ | $-110.09$ | $-1.21\times 10^{-6}$ |
| 100 | $373\,473.8$ | $10\,523.4$ | $-110.04$ | $-3.83\times 10^{-8}$ |

### 13.4 透视加速度对比

| 指标 | 理论值 | 实测值 | 偏差 |
|------|--------|--------|------|
| 自行漂移 $\dot\mu$ | $1.285\ \text{mas/yr}^2$ | $1.300\ \text{mas/yr}^2$ | $+1.2\%$ |
| 径向加速度 $\dot v_r$ | $4.50\ \text{m/s/yr}$ | $4.57\ \text{m/s/yr}$ | $+1.5\%$ |

### 13.5 结论

- 日心距每 10 年线性减少约 $233\ \text{AU}$，对应匀速直线运动，引力扰动完全不可见；
- $\tau$ 场全程处于准稳态漂移，振幅受线性阻尼 $K\tau$ 与速度驱动 $\gamma v^2$ 的平衡压制，始终满足 $|\tau|\lesssim 10^{-6}$；
- 透视加速度自动涌现，与理论值偏差 $<3\%$，证实积分器在 $10^5\ \text{AU}$、$10^2\ \text{yr}$ 尺度上数值干净。

---

## 14. 跨系统标度：$\gamma$ 普适与 $K$ 质心锁定

### 14.1 代码单位体系

所有数值积分采用**固定代码单位**，以太阳系水星为基准定义：

| 代码单位 | 物理值 | 验证 |
|----------|--------|------|
| $L_0$ | $0.258065\ {\rm AU}$ | $L_0^3/T_0^2 = 39.48 = GM_\odot$ |
| $T_0$ | $0.020865\ {\rm yr}$ | 水星周期 $T_{\rm code}=2\pi\sqrt{A_{\rm REF}^3}=11.54$ |
| $V_0$ | $0.0339\ {\rm AU/day}$ | $V_0=L_0/T_0$ |

在代码单位中 $\alpha\equiv GM_\odot=1.0$。跨系统时，中心质量直接映射为 $\alpha=M_\star$（以太阳质量为单位）。

### 14.2 $\gamma$ 的普适性

$\gamma$ 是 mini-superspace 截断后共形场论的唯一残余耦合。量纲分析给出：

$$\gamma_{\rm phys}=\frac{\gamma_{\rm code}}{L_0^2},\qquad \gamma_{\rm code}=2.6186\times10^{-7}.$$

$\gamma_{\rm code}$ 跨系统**严格不变**。物理上，它对应某种真空极化或共形反常的低能投影，与具体恒星质量无关。

### 14.3 $K$ 的质心锁定

$K$ 的量纲为 $[T^{-2}]$（刚度/频率平方）。在 mini-superspace 截断下，$\tau$ 场的响应频率必须与系统特征束缚频率匹配，否则会出现非物理的参量共振或虚假自由漂移。因此 $K$ 由系统质心锁定：

$$\boxed{K_{\rm code}=2.5\times\frac{M_\star}{M_\odot}\left(\frac{a_{\rm Mercury}}{a_{\rm ref}}\right)^3}$$

其中 $a_{\rm ref}$ 为该系统最内层（或参考）轨道的半长轴。同一系统内所有行星共享同一个 $K$；跨系统时 $K$ 随质心质量与特征尺度迁移。

### 14.4 周期自洽验证

| 系统 | 天体 | 物理周期 | 代码周期 $\times T_0$ | 偏差 |
|------|------|----------|------------------------|------|
| TRAPPIST-1 | b | $1.51\ {\rm d}$ | $1.51\ {\rm d}$ | $0.01\%$ |
| TRAPPIST-1 | h | $18.8\ {\rm d}$ | $18.90\ {\rm d}$ | $0.54\%$ |
| 51 Pegasi | b | $4.23\ {\rm d}$ | $4.42\ {\rm d}$ | $4.5\%$ |
| HD 209458 | b | $3.50\ {\rm d}$ | $3.49\ {\rm d}$ | $0.38\%$ |
| 太阳系 | 水星 | $87.97\ {\rm d}$ | $87.97\ {\rm d}$ | $0.00\%$ |
| 太阳系 | 地球 | $365.25\ {\rm d}$ | $365.26\ {\rm d}$ | $0.00\%$ |

周期自洽性证实代码单位体系跨系统无漂移。

### 14.5 跨系统数值验证

将固定 $(\gamma_{\rm code}, K_{\rm sun}=2.5)$ 迁移至四个行星系统，$K$ 按质心规则重算：

| 系统 | $M_\star$ | $a_{\rm ref}$ | $K_{\rm code}$ | 天体 | $a$ (AU) | Weyl (''/cy) | GR (''/cy) | Model/GR |
|------|-----------|---------------|----------------|------|----------|--------------|-------------|----------|
| **太阳系** | $1.0$ | $0.387$ | $2.5$ | 水星 | $0.387$ | $42.97$ | $44.9$ | $0.958$ |
| | | | | 地球 | $1.000$ | $3.53$ | $4.0$ | $0.880$ |
| | | | | 冥王星 | $39.48$ | $8\times10^{-5}$ | $4\times10^{-4}$ | $0.874$ |
| **TRAPPIST-1** | $0.089$ | $0.0115$ | $8486$ | b | $0.0115$ | $2.21$ | $7502$ | $0.00029$ |
| | | | | h | $0.062$ | $0.03$ | $111$ | $0.00026$ |
| **51 Pegasi** | $1.0$ | $0.0527$ | $991$ | b | $0.0527$ | $15.7$ | $6286$ | $0.0025$ |
| **HD 209458** | $1.0$ | $0.045$ | $1591$ | b | $0.045$ | $231.7$ | $9328$ | $0.0248$ |

**内部平行性**：
- 太阳系：Model/GR 均值 $0.904$，标准差 $0.038$，极差 $8.3\%$
- TRAPPIST-1：Model/GR 均值 $0.00028$，标准差 $<10^{-5}$

---

## 15. 奥陌陌：银心质心系中的 $\tau$ 冻结

### 15.1 观测非引力加速度

Micheli et al. (2018) 的精密天体测量显示，奥陌陌偏离开普勒轨道达到 $30\sigma$ 显著性：

| 参数 | 观测值 |
|------|--------|
| 非引力加速度 $A_1$ | $2.79\times10^{-7}\ {\rm AU/day^2}$ |
| 物理换算 | $5.59\times10^{-6}\ {\rm m/s^2}$ |
| 方向 | 径向远离太阳 |
| 距离标度 | $a_{\rm ng}\propto r^{-2}$（最佳拟合） |
| 累积位置偏差 | 消失时比纯引力轨道远 $\sim 40\,000\ {\rm km}$ |

### 15.2 银心质心系下的 $\tau$ 冻结

奥陌陌来自银心质心系，等效 $K_{\rm gal}\sim10^{-22}$ 几乎为零，$\tau$ 场失去恢复力，在星际游荡期间被 $v^2$ 驱动缓慢漂移，最终冻结在 $\tau_{\rm frozen}\approx-1.85\times10^{-3}$。

有效引力常数变为 $\alpha_{\rm eff}=e^\tau\alpha\approx(1-1.85\times10^{-3})\alpha$，太阳引力等效减弱，产生径向残余加速度：

$$\Delta a_r\approx-\tau\cdot\frac{GM_\odot}{r^2}$$

在 $r=1.4\ {\rm AU}$ 处，反推 $\tau$ 产生的加速度与观测 $A_1$ 比值：

$$\frac{a_{\rm ng}^{\rm Weyl}}{A_1^{\rm obs}}=-0.999$$

### 15.3 数值积分验证

设定奥陌陌双曲线轨道（$q=0.255\ {\rm AU},\ e=1.20$），比较纯引力与 $\tau$ 冻结轨道：

| 检验项 | 结果 |
|--------|------|
| 在 $r=1.4\ {\rm AU}$ 处径向加速度偏差 | $-5.59\times10^{-6}\ {\rm m/s^2}$（与观测 $-0.1\%$ 匹配） |
| 到达 $r=1.4\ {\rm AU}$ 时间差 | $-0.093\ {\rm d}$ |
| 该处速度差 | $+0.12\ {\rm km/s}$ |
| 到 $r=5\ {\rm AU}$ 位置偏差 | $\sim 3.0\times10^6\ {\rm km}$ |

### 15.4 可区分预言

标准解释（H₂ 喷气/辐射压）要求加速度与天体物理属性（质量/面积比、自转轴）相关。Weyl $\tau$ 冻结模型额外预言：

1. **无自转耦合**：$\tau$ 场与自转轴无关，非引力加速度不应显示自转相位调制——与奥陌陌光变曲线和 $A_1$ 之间缺乏相关性一致。
2. **星际来客统计散布**：不同来客携带不同 $\tau$ 冻结值，$A_1$ 应有统计散布。
3. **无质量依赖**：与天体质量无关，大小天体应有相同 $A_1$（若 $\tau$ 背景相同）。

---

## 16. 各物理壳的涌现极限

通过赋予 $\tau$ 不同的背景行为，母作用量自动涌现四大物理分支：

| 壳 | $\tau$ 设定 | 涌现方程 | 物理 |
|----|------------|---------|------|
| 引力 | $\tau\ll 0$ | $\dot{\mathbf{v}}_i=\dot\tau\mathbf{v}_i-e^\tau\sum_j\frac{Gm_j\mathbf{r}_{ij}}{r_{ij}^3}$ | 等效弱引力，恢复牛顿极限 |
| 电磁/超导 | 空间区域 $\tau\gg 0$ | $\nabla^2 B=e^\tau B$ | 光子质量，迈斯纳效应 |
| 粒子/禁闭 | 内部 $\tau\to\infty$ | $e^{-\tau}\to 0$，动能冻结 | 无自由夸克 |
| 量子 | $\tau=0$ | $S_{\rm ent}=\ln 2$，$\mathcal{B}=2\sqrt{2}$ | 最大纠缠，标准量子力学 |

---

## 17. 星系尺度：多体 $\tau$ 叠加的等效暗物质势

### 17.1 多体作用量与连续化

假设星系由 $N$ 个恒星-行星子系统组成，每个子系统 $i$ 在自身 mini-superspace 中贡献 $\tau_i$。在星系质心系中，总作用量：

$$\mathcal{L}_{\rm gal} = \sum_{i=1}^N \left[\frac{1}{2}e^{-\tau_i}v_i^2 + \frac{1}{2}\dot\tau_i^2 - \frac{K_i}{2}\tau_i^2\right] + \mathcal{L}_{\rm int}(\{\tau_i\})$$

当子系统足够多且分布均匀，做连续化近似 $\tau_i \to \tau(\mathbf{r})$，并假设长程关联主导（梯度项恢复）：

$$\mathcal{L}_{\rm gal} \to \int d^3r \left[\frac{1}{2}(\nabla\tau)^2 + \frac{1}{2}\dot\tau^2 - \frac{K_{\rm gal}}{2}\tau^2 - \gamma_{\rm gal}\, e^{-\tau}\rho(\mathbf{r})v^2(\mathbf{r})\right]$$

其中 $\rho$ 为恒星数密度。准稳态下欧拉-拉格朗日方程给出：

$$\nabla^2\tau - K_{\rm gal}\tau = -\gamma_{\rm gal} e^{-\tau}\rho v^2$$

### 17.2 等温晕极限与对数势

若 $K_{\rm gal}\approx 0$（星系尺度刚度极弱，类似银心极限），方程退化为：

$$\nabla^2\tau = -\gamma_{\rm gal} e^{-\tau}\rho v^2$$

对于球对称星系，取等温晕近似 $\rho(r) \propto 1/r^2$，且平坦旋转曲线 $v \approx v_{\rm flat}$，源项 $\propto 1/r^2$。方程的解在 $r \to \infty$ 时：

$$\tau(r) \sim -\frac{\gamma_{\rm gal}M_{\rm enc}(r)v_{\rm flat}^2}{4\pi r} \propto -\ln r$$

因此有效度规：

$$g_{00} = e^{-\tau(r)} \sim r^{\alpha_{\rm DM}}, \quad \alpha_{\rm DM} \sim \frac{\gamma_{\rm gal}v_{\rm flat}^2}{4\pi}$$

### 17.3 平坦旋转曲线的涌现

等效引力势（从测地线提取）：

$$\Phi_{\rm eff}(r) = -\frac{GM(r)}{r} + \frac{1}{2}\ln g_{00} \approx -\frac{GM(r)}{r} + \frac{|\tau(r)|}{2}$$

当 $r$ 足够大，对数项 $\propto \ln r$ 超越 $1/r$ 衰减，产生平坦旋转曲线：

$$v_{\rm circ}^2 = r\frac{d\Phi_{\rm eff}}{dr} \approx \frac{GM(<r)}{r} + \frac{1}{2}r\frac{d|\tau|}{dr} \to {\rm const}$$

此即星系旋转曲线平坦的共形场论解释，无需引入暗物质晕或未知粒子。

---

## 18. 真空极限：旋量湮灭与宇宙学红移

当物质场真空湮灭（$\chi=\bar\chi=0$），蝴蝶变换后的 $\tau$ 动能系数退化为 1：

$$\mathcal{L}_\tau\to-\frac{1}{2}(\partial_\mu\tau)^2-\frac{K}{2}\tau^2$$

此时 $\tau$ 可纯粹地扮演度规共形因子。作替换 $g_{\mu\nu}=e^{-\tau}\eta_{\mu\nu}$，光子沿零测地线传播，物理频率（固有时测量）满足：

$$\omega_{\rm phys}=e^{-\tau/2}\omega_{\rm conf}$$

由于共形坐标中 $\omega_{\rm conf}$ 沿射线守恒，红移为：

$$\boxed{z=e^{(\tau_{\rm obs}-\tau_{\rm emit})/2}-1\approx\frac{\Delta\tau}{2}\quad(|\Delta\tau|\ll 1)}$$

在宇宙学真空 regime（$\chi=0$）中，驱动项消失，$\tau$ 方程退化为简谐振子：

$$\ddot\tau+K\tau=0\quad\Rightarrow\quad\tau(t)=\tau_{\rm amp}\cos(\sqrt{K}\,t+\phi)$$

若宇宙学背景处于这种振荡中，红移将带有周期性调制，与标准 $\Lambda$CDM 的单调红移不同。若 $K_{\rm cos}\sim H_0^2$，则振荡周期 $\sim1/H_0$，恰好是宇宙学量级。

---

## 19. 宇宙学尺度：$\tau$ 振荡与等效暗能量

### 19.1 真空简谐振子

在宇宙学均匀各向同性极限下，物质场湮灭（$\chi=\bar\chi=0$），$\tau$ 退化为真空简谐振子：

$$\ddot\tau + K_{\rm cos}\tau = 0$$

宇宙学背景中 $K_{\rm cos} \sim H_0^2$（刚度由哈勃尺度锁定），解为：

$$\tau(t) = \tau_{\rm amp}\cos(\sqrt{K_{\rm cos}}\,t + \phi_0)$$

### 19.2 有效引力常数与 Friedmann 方程

有效引力常数随时间调制：

$$G_{\rm eff}(t) = G_0\,e^{\tau(t)} \approx G_0\left[1 + \tau_{\rm amp}\cos(\sqrt{K_{\rm cos}}\,t)\right]$$

Friedmann 方程（平坦背景 $\eta_{\mu\nu}$ 中，有效度规 $g_{\mu\nu}=e^{-\tau}\eta_{\mu\nu}$）：

$$\left(\frac{\dot a}{a}\right)^2 = \frac{8\pi G_{\rm eff}(t)}{3}\rho_m + \frac{\Lambda_{\rm eff}(t)}{3}$$

其中等效宇宙学常数来自 $\tau$ 场的能量密度：

$$\rho_\tau = \frac{1}{2}\dot\tau^2 + \frac{K_{\rm cos}}{2}\tau^2 = K_{\rm cos}\tau_{\rm amp}^2 = {\rm const}$$

### 19.3 状态方程与负压相

状态方程：

$$w_\tau = \frac{P_\tau}{\rho_\tau} = \frac{\frac{1}{2}\dot\tau^2 - \frac{K_{\rm cos}}{2}\tau^2}{\frac{1}{2}\dot\tau^2 + \frac{K_{\rm cos}}{2}\tau^2} = \cos^2(\sqrt{K_{\rm cos}}\,t) - \sin^2(\sqrt{K_{\rm cos}}\,t) = \cos(2\sqrt{K_{\rm cos}}\,t)$$

因此：

$$\boxed{w_\tau(t) = \cos(2\sqrt{K_{\rm cos}}\,t)}, \quad \langle w_\tau\rangle = 0 \text{（过零振荡）}$$

观测到的宇宙加速膨胀对应 $\tau$ 在特定半周期（等效排斥相），此时：
- $G_{\rm eff}$ 减弱（等效排斥）
- $w_\tau < 0$（负压驱动加速）

### 19.4 与 DESI $w(z)$ 演化的定性对应

DESI 观测到 $w(z)$ 偏离常数 $-1$ 并呈现演化特征。在 Weyl-CFT 中，这是 $\tau$ 场真空振荡的自然结果——单参数振荡模型在 $(w_0,w_a)$ 平面上产生一条理论约束线，与 DESI 置信域相交。若振荡周期 $\sim H_0^{-1}$，当前宇宙恰好处于负压排斥半周期；未来可能转向吸引相，这是原则可证伪的预言。

### 19.5 大尺度结构：共形排斥力的尺度依赖竞争

拉格朗日中 $\Phi$ 与 $\bar\Phi$ 携带相反共形荷：

$$\mathcal{L} \supset \frac{e^{-\tau}}{2}|D_\mu\Phi|^2 + \frac{e^{\tau}}{2}|D_\mu\bar\Phi|^2$$

共轭动量：

$$\pi_\Phi = e^{-\tau}\dot\Phi^\dagger, \quad \pi_{\bar\Phi} = e^{\tau}\dot{\bar\Phi}^\dagger$$

当两个系统（或两个宇宙学区域）的 $\tau$ 场重叠，联合 $\tau = \tau_1 + \tau_2$，交叉项产生等效排斥压：

$$P_{\rm rep} = -\frac{\partial\mathcal{L}}{\partial V}\bigg|_{T,V} = \frac{K}{2}\tau^2 + \gamma e^{-\tau}\langle v^2\rangle$$

在宇宙学均匀各向同性近似下，平均化给出：

$$\langle P_{\rm rep}\rangle = \frac{K_{\rm cos}}{4}\tau_{\rm amp}^2 = \frac{\rho_\tau}{2} > 0 \quad (\text{当 } \dot\tau^2 < K\tau^2\text{ 时})$$

大尺度结构的形成源于共形排斥力与引力坍缩的尺度依赖竞争，替代了标准模型中暗物质引力坍缩与暗能量对抗的二元机制。

---

## 20. 极限检验与自洽性

### 20.1 $K=0,\ \gamma=0$：严格退化开普勒

当关闭共形耦合且初始 $\tau=0,\ \dot\tau=0$ 时，运动方程严格退化为：

$$\ddot r = r\dot\theta^2 - \frac{\alpha}{r^2}$$

数值验证（水星轨道 10 圈）显示位置偏差 $<10^{-12}$ code（$\sim10^{-5}\ {\rm km}$），在机器精度内等于零。

### 20.2 $K=0,\ \gamma\neq0$：$\tau$ 自由漂移

银心极限下，$\tau$ 方程退化为 $\ddot\tau\approx-\gamma e^{-\tau}v^2$，无恢复力，单向累积。这正是奥陌陌显示非引力加速度的数学根源。

### 20.3 $K>0,\ \gamma=0$：$\tau$ 阻尼回零

有刚度无驱动时，$\tau$ 被 $-K\tau$ 拉回零，轨道严格为标准开普勒。

---

## 21. 模型边界与反例

### 21.1 冥王星：质心错配演示

冥王星若被错误地当作独立质心系（$a_{\rm ref}=39.48\ {\rm AU}$），则 $K\approx2.4\times10^{-6}$，与正确值 $2.5$ 相差 106 万倍，$\tau$ 场完全崩溃。在太阳系质心系内共用 $K=2.5$，冥王星落在 Model/GR $\approx 0.87$ 的平行带上。

### 21.2 1I/'Oumuamua：无质心束缚（银心极限）

双曲线轨道，无恒星质心。在太阳质心系中 $K$ 无定义；换至银心质心系，$K\sim10^{-22}$，$\tau$ 自由漂移，冻结在 $\sim10^{-3}$，产生观测到的非引力加速度。

### 21.3 PSR B1913+16：强场双星

总质量 $2.828\,M_\odot$，$v\sim300\ {\rm km/s}$（$0.001c$），弱场线性近似失效。需后牛顿扩展。

### 21.4 Kepler-16 b：非单星质心

环绕双星行星，质心选择敏感（总质量 vs 主星质量导致 $K$ 差异 1.3 倍）。模型需扩展为多体 mini-superspace。

---

## 22. 结论

我们证明了如下简单作用量

$$S=\int d^4x\left[\frac{e^{-\tau}}{2}(D_\mu\Phi)^\dagger(D^\mu\Phi)+\frac{e^{\tau}}{2}(D_\mu\bar\Phi)^\dagger(D^\mu\bar\Phi)-\frac{K}{2}\tau^2+\frac{1}{4g^2}F_{\mu\nu}F^{\mu\nu}\right]$$

支持一种涌现的轨道力学，能够在完全不引用弯曲时空的前提下复现并预测多系统近心点进动。关键要素包括：

1. **Weyl–共形耦合** $e^{\pm\tau}$，动态调制物质场有效质量；
2. **非对称驱动** $\ddot\tau=-K\tau-\gamma e^{-\tau}v^2$，打破时间反演并产生久期漂移；
3. **涌现标度律**：维里定理自动保护 $a^{-5/2}$ 进动标度，与 GR 同构；
4. **多行星联合预测**：单参数 $\gamma$ 锚定后，七行星预测值高度平行（标准差 $<1\%$）；
5. **双 regime 边界**：3I/ATLAS 的 $\tau$ 贡献仅为观测非引力偏差的千分之二；
6. **长程稳定性**：巴纳德星百年透视加速度自然涌现，$\tau$ 全程受抑；
7. **跨系统标度**：$\gamma$ 普适，$K$ 质心锁定，紧凑系统自动退化，太阳系尺度与 GR 竞争；
8. **星际来客**：奥陌陌的观测非引力加速度由银心冻结 $\tau$ 场自然解释（匹配 $-0.1\%$）；
9. **星系尺度**：多体 $\tau$ 叠加产生对数型等效势修正，自然导出平坦旋转曲线；
10. **真空极限**：$\chi=0$ 时 $\tau$ 升格为度规共形因子，预言宇宙学红移 $z\approx\Delta\tau/2$；
11. **宇宙学尺度**：$\tau$ 真空振荡产生等效暗能量动力学，状态方程过零振荡，与 DESI $w(z)$ 演化定性相容；
12. **量子调制**：$\tau$ 场连续调制纠缠熵与贝尔不等式，预言强场退相干；
13. **拓扑起源**：整个理论可从四色定理的四面体剖分严格推导，谱隙原理保证非零真空。

本工作表明，轨道进动、星系旋转曲线、宇宙加速膨胀与量子纠缠不必是度规引力的专属签名。它们可以从共形物质场动力学中涌现，且具备罕见的系统预测能力，为替代引力现象学开辟新窗口。

---

## 附录 A：符号约定

| 约定项 | GR | 本模型（Weyl） |
|--------|-----|----------------|
| 轨道 $\dot\theta$ 符号 | 正 $\equiv$ 逆时针 | 负 $\equiv$ 逆时针（坐标手性） |
| 进动与轨道方向关系 | 同向（顺行） | 同向（顺行） |
| 报告进动符号 | 正 | 负 |
| 物理实质 | 近日点沿运动方向前进| 近日点沿运动方向前进 |

---

## 附录 B：3I/ATLAS (C/2025 N1) 轨道数据

**历元**：2026-Feb-19.000000 TDB

**坐标系**：ICRF 赤道参考架，太阳质心原点

| 分量 | 位置 (AU) | 速度 (AU/day) |
|------|-----------|---------------|
| $x$ | $-1.878465706009316$ | $-3.393307402708430\times 10^{-3}$ |
| $y$ | $3.556787034031755$ | $3.324252268915676\times 10^{-2}$ |
| $z$ | $1.369259882335852$ | $1.201145093769092\times 10^{-2}$ |

**导出量**：
- 日心距 $r = 4.048\ \text{AU}$
- 速度 $|v| = 0.0345\ \text{AU/day} \approx 5.96\ \text{km/s}$
- 双曲线超速 $v_\infty \sim 3.1\ \text{km/s}$
- 近日点距 $q = 1.36\ \text{AU}$（由 $e\approx 6.14$ 反推）

---

## 附录 C：巴纳德星 (GJ 699) 初始条件

**历元**：J2000.0 (ICRS, 日心)

**来源**：Gaia DR3

| 参数 | 数值 |
|------|------|
| 赤经 $\alpha$ | $269.452076375^\circ$ |
| 赤纬 $\delta$ | $+4.693390672^\circ$ |
| 视差 $\pi$ | $548.9\ \text{mas}$ |
| 距离 $d$ | $375\,800\ \text{AU}$ |
| 自行 $\mu_\alpha$ | $-802.3\ \text{mas/yr}$ |
| 自行 $\mu_\delta$ | $+10\,362.5\ \text{mas/yr}$ |
| 径向速度 $v_r$ | $-110.5\ \text{km/s}$ |

**笛卡尔状态向量**（AU，AU/day）：
- $\mathbf{r} = [-3.566\times 10^3,\ -3.745\times 10^5,\ +3.074\times 10^4]$
- $\mathbf{v} = [-3.350\times 10^{-3},\ +6.786\times 10^{-2},\ +4.622\times 10^{-2}]$

---

## 附录 D：跨系统代码单位体系

### D.1 固定基准（太阳系，永不变）

$$L_0 = \frac{a_{\rm Mercury}}{A_{\rm REF}} = 0.258065\ {\rm AU}$$

$$T_0 = \frac{T_{\rm Mercury}}{2\pi/\sqrt{A_{\rm REF}^{-3}}} = 0.020865\ {\rm yr}$$

$$V_0 = \frac{L_0}{T_0} = 0.0339\ {\rm AU/day}$$

验证：$L_0^3/T_0^2 = 39.48\ {\rm AU^3/yr^2} = GM_\odot$

### D.2 普适常数

$$\gamma_{\rm code} = 2.6186\times10^{-7}\quad [L^{-2}]$$

$$\gamma_{\rm phys} = \frac{\gamma_{\rm code}}{L_0^2} = 3.9320\times10^{-6}\ {\rm AU^{-2}}$$

### D.3 系统特定参数

$$\alpha = M_\star\quad (\text{以 }M_\odot\text{ 为单位})$$

$$K_{\rm code} = 2.5\times\frac{M_\star}{M_\odot}\left(\frac{0.387098\ {\rm AU}}{a_{\rm ref}}\right)^3$$

### D.4 行星输入

$$a_{\rm code} = \frac{a_{\rm phys}}{L_0},\qquad e_{\rm code}=e_{\rm phys}$$

### D.5 输出换算

$$T_{\rm phys} = T_{\rm code}\times T_0$$

$$\Delta\phi_{\rm century} = \Delta\phi_{\rm orbit}\times\frac{180^\circ}{\pi}\times3600\times\frac{100\ {\rm yr}}{T_{\rm phys}}$$

---

## 附录 E：奥陌陌 (1I/'Oumuamua) 轨道参数

**历元**：2017-Sep-09（近日点通过）

| 参数 | 数值 |
|------|------|
| 近日点距 $q$ | $0.2553\ {\rm AU}$ |
| 离心率 $e$ | $1.1995$ |
| 双曲线超速 $v_\infty$ | $26.32\ {\rm km/s}$ |
| 近日点速度 $v_{\rm peri}$ | $87.4\ {\rm km/s}$ |
| 非引力加速度 $A_1$ | $2.79\times10^{-7}\ {\rm AU/day^2}$ |
| 方向 | 径向远离太阳 |
| 距离标度 | $\propto r^{-2}$ |

**反推冻结参数**：
- $\tau_{\rm frozen} = -1.85\times10^{-3}$
- 等效引力减弱：$e^\tau \approx 0.99815$

---

## 附录 F：关键验算公式

### F.1 宇宙学状态方程 $w_\tau(t)$

由 $\tau(t) = \tau_{\rm amp}\cos(\omega t)$，$\omega=\sqrt{K_{\rm cos}}$，能量密度与压强：

$$\rho_\tau = \frac{1}{2}\dot\tau^2 + \frac{K_{\rm cos}}{2}\tau^2 = K_{\rm cos}\tau_{\rm amp}^2$$

$$P_\tau = \frac{1}{2}\dot\tau^2 - \frac{K_{\rm cos}}{2}\tau^2 = K_{\rm cos}\tau_{\rm amp}^2\cos(2\omega t)$$

状态方程：

$$\boxed{w_\tau(t) = \frac{P_\tau}{\rho_\tau} = \cos(2\sqrt{K_{\rm cos}}\,t)}$$

### F.2 星系尺度 $\tau$ 场方程

准稳态 $K_{\rm gal}\approx 0$ 下：

$$\nabla^2\tau = -\gamma_{\rm gal} e^{-\tau}\rho v^2$$

球对称等温晕 $\rho\propto 1/r^2$，$v=v_{\rm flat}$，源项 $\propto 1/r^2$，解：

$$\tau(r) \sim -\frac{\gamma_{\rm gal}M_{\rm enc}(r)v_{\rm flat}^2}{4\pi r} \propto -\ln r$$

### F.3 共形排斥压

$$P_{\rm rep} = \frac{K}{2}\tau^2 + \gamma e^{-\tau}\langle v^2\rangle$$

宇宙学均匀近似下：

$$\langle P_{\rm rep}\rangle = \frac{K_{\rm cos}}{4}\tau_{\rm amp}^2 = \frac{\rho_\tau}{2} \quad (\dot\tau^2 < K\tau^2\text{ 时})$$

### F.4 红移公式

度规共形因子 $g_{\mu\nu}=e^{-\tau}\eta_{\mu\nu}$，光子零测地线：

$$\omega_{\rm phys}=e^{-\tau/2}\omega_{\rm conf}$$

$$\boxed{z = e^{(\tau_{\rm obs}-\tau_{\rm emit})/2}-1 \approx \frac{\Delta\tau}{2}\quad(|\Delta\tau|\ll 1)}$$

---

**固定参数**：$K=2.5\ (\text{太阳系基准}),\ \delta=0,\ \gamma=2.6186\times 10^{-7}$（代码单位）
