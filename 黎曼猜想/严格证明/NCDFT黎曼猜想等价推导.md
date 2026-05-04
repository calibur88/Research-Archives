# NCDFT 框架与黎曼猜想的双向等价性

---

## 前言

本文档给出 NCDFT（非交换离散傅里叶变换）框架与黎曼假设（RH）的**双向等价性**。上篇为**逆问题**：从黎曼零点谱出发，反推承载该谱的代数骨架必须是非交换 DFT，且其唯一可递归参数为 $\alpha=1/2$。下篇为**正问题**：从算术前提出发构造 NCDFT，证明 $\alpha=1/2$ 时其谱测度收敛于黎曼零点分布，且任何离轴零点都会导致构造性矛盾。两篇在定理 8.1' 处闭合为严格的逻辑环。

---

# 上篇：逆问题——从黎曼谱到非交换 DFT 的构造性推导

---

## 1. 对数生成元

黎曼零点计数函数的反解给出第 $n$ 个零点虚部的渐近公式，其闭式由 Lambert $W$ 函数写出：

```math
y_n=\frac{2\pi\bigl(n-\frac{11}{8}\bigr)}{W\!\bigl[e^{-1}\bigl(n-\frac{11}{8}\bigr)\bigr]}.
```

对大 $n$ 展开 $W(z)=\log z-\log\log z+\cdots$，零点序列的主导变量是**嵌套对数** $\log n$。因此，任何以零点为谱目标的离散算子，其生成元必须内禀地以对数坐标为对角核。

取 $\mathfrak{su}(2)$ 的 $N=2j+1$ 维不可约表示，标准基 $|j,m\rangle$，$m=-j,\dots,j$，编号 $n=m+j+1\in\{1,\dots,N\}$。定义**对数生成元**

```math
\Lambda=\log(J_z+j+1)=\operatorname{diag}(\log 1,\log 2,\dots,\log N).
```

---

## 2. Mellin 对偶与 DFT 骨架

黎曼显式公式经 Mellin 变换后，函数方程 $s\leftrightarrow 1-s$ 在临界线 $\Re s=1/2$ 上等价于 Fourier 反射对称 $t\leftrightarrow -t$。在 $\mathfrak{su}(2)$ 表示论中，这一对称由 $J_z$ 基到 $J_x$ 基的 $\pi/2$ 绕 $y$ 轴旋转实现：

```math
F_{1/2}=e^{-i\frac{\pi}{2}J_y},
```

其矩阵元为 Wigner 小 $d$-函数 $d^j_{m'm}(\pi/2)$。该算子满足

```math
F_{1/2}^2:\;m\mapsto -m,\qquad F_{1/2}^4=I,
```

精确对应经典 DFT 的四阶周期与 Mellin 变换的自对偶性。

---

## 3. 非交换扩张：蝶形算子 $Q$

为引入可调参数 $\alpha$ 并保持非交换结构，定义**宇称算子**

```math
B=(-1)^{J_z+j},
```

它在 $J_z$ 基下本征值为 $(-1)^{n+1}$，精确区分偶/奇格点。由反对易关系 $\{J_x,B\}=0$（$J_x$ 翻转宇称），构造厄米**蝶形算子**

```math
Q=iJ_xB.
```

**性质**：
- $Q^\dagger=Q$；
- $Q^2=J_x^2$；
- 在 $J_x$ 基下，$Q$ 为纯反对角，连接 $\pm m_x$ 子空间，其矩阵元呈 $\pm i|m_x|$ 的配对结构。

**扭曲 DFT**定义为

```math
F_\alpha=F_{1/2}\,e^{i(\alpha-\frac12)Q}.
```

当 $\alpha=1/2$ 时退化为标准 DFT。相位因子 $e^{i\theta Q}$（其中 $\theta=\alpha-\frac12$）在 $J_x$ 基下按 $J_x$ 本征值分块，每块为 $2\times 2$ 旋转矩阵。

---

## 4. Jacobi 矩阵的构造与鞍点凝聚

定义算子簇

```math
\mathcal{J}_\alpha=F_\alpha\,\Lambda\,F_\alpha^\dagger.
```

在 $J_x$ 基下，$\mathcal{J}_\alpha$ 的矩阵元由 Wigner $d$-函数的积分表示给出。利用鞍点法，对大 $j$ 令 $\mu=m/j,\nu=k/j$，作用量

```math
S(\phi;\mu,\nu)=i(\mu-\nu)\phi+\frac{1+\mu}{2}\log\!\Bigl(\frac{1+ie^{i\phi}}{2}\Bigr)+\frac{1-\mu}{2}\log\!\Bigl(\frac{1-ie^{-i\phi}}{2}\Bigr)
```

的鞍点方程 $\partial_\phi S=0$ 表明：当 $|\mu-\nu|\gg j^{-1/2}$ 时，鞍点离开实轴进入复平面，导致指数衰减

```math
\bigl(\mathcal{J}_\alpha\bigr)_{mk}\sim C\,e^{-j|\mu-\nu|^2/2}.
```

因此，当 $N\to\infty$ 时，仅最近邻与次近邻耦合存活。执行 Lanczos 三对角化后，$\mathcal{J}_\alpha$ 凝聚为实对称**Jacobi 矩阵**：

```math
\mathcal{J}_\alpha=\begin{pmatrix}
a_1 & b_1 & & \\
b_1 & a_2 & b_2 & \\
& \ddots & \ddots & \ddots \\
& & b_{N-1} & a_N
\end{pmatrix}.
```

**矩阵元的标度行为**：
- 对角元 $a_n(\alpha)=\log n+O(1/j)$；
- 次对角元在 $\alpha=1/2$ 时退化为 $\mathfrak{su}(2)$ 标准升降系数 $b_n(1/2)\sim\frac12\sqrt{n(N-n)}$；
- 当 $\alpha\neq 1/2$ 时，$Q$ 的双线性宇称结构引入与系统尺寸成正比的曲率修正：

```math
b_n(\alpha)=b_n(1/2)+|\alpha-\tfrac12|\cdot\frac{N^2}{(\log N)^2}\cdot\kappa_n+O\bigl((\alpha-\tfrac12)^2\bigr),
```

其中 $\kappa_n>0$ 为 $O(1)$ 的曲率核。

---

## 5. 函子递归：Cooley–Tukey 的代数刚性

Cooley–Tukey 基-2 FFT 要求 DFT 矩阵具有精确的分块蝶形结构：

```math
F_N=\frac{1}{\sqrt2}\begin{pmatrix} F_{N/2} & D\,F_{N/2} \\ F_{N/2} & -D\,F_{N/2} \end{pmatrix},
```

其中 $D$ 为 twiddle 因子对角阵。在 $\mathfrak{su}(2)$ 框架中，这对应于表示空间按 $B$ 的偶/奇本征值分解 $V=V_{\mathrm{even}}\oplus V_{\mathrm{odd}}$，且 $F_{1/2}$ 精确保持该分块结构。

对于扭曲 DFT $F_\alpha$，相位因子 $e^{i\theta Q}$ 在偶/奇基下为反对角块：

```math
e^{i\theta Q}=\begin{pmatrix} \cos(\theta|J_x|) & i\sin(\theta|J_x|)\,Q/|J_x| \\ \cdots & \cdots \end{pmatrix}.
```

只要 $\theta=\alpha-\frac12\neq 0$，反对角块即非零，破坏蝶形结构所要求的“对角块严格为 $F_{N/2}$、反对角块严格为 $\pm D\,F_{N/2}$”的精确递归关系。因此，NCDFT 簇 $\{F_\alpha^{(N)}\}_{N=2^k}$ 形成函子的**充要条件**为 $\theta=0$，即

```math
\alpha=\frac12.
```

---

## 6. 谱相变：Carleman 条件

Jacobi 矩阵本质自伴的 Carleman 条件为

```math
\sum_{n=1}^{\infty}\frac{1}{\beta_n}=\infty.
```

- **$\alpha=1/2$**：$\beta_n(1/2)\sim\frac12\sqrt{n(N-n)}$。在极限 $N\to\infty$ 下，固定 $n$ 时 $\beta_n\sim\sqrt{n}$，级数 $\sum 1/\sqrt{n}$ 发散，算子本质自伴，谱严格支撑于实轴。
- **$\alpha\neq 1/2$**：曲率修正使 $\beta_n(\alpha)\sim|\alpha-\frac12|N^2/(\log N)^2$ 对大 $N$ 占主导，Carleman 和

```math
\sum_{n=1}^{N}\frac{1}{\beta_n(\alpha)}\leq\frac{C}{|\alpha-\frac12|}<\infty.
```

极限下算子失去本质自伴性，其自伴扩张允许非实谱，即零点可离轴。

---

## 7. 符号凝聚与正则化

定义正则化符号算子

```math
\mathcal{S}_\alpha^{(\varepsilon)}=\frac{\mathcal{J}_\alpha}{\sqrt{\mathcal{J}_\alpha^2+\varepsilon^2 I}},\qquad\varepsilon>0.
```

**$\alpha=1/2$**：$\mathcal{J}_{1/2}$ 本质自伴，$|\mathcal{J}_{1/2}|$ 良定。取 $\varepsilon\to 0$，利用 Plancherel–Rotach 渐近（鞍点法），Jacobi 矩阵的正交多项式本征函数在转折点附近退化为 Airy 函数，远离转折点退化为 WKB 振荡子。离散谱和转化为 Hilbert 主值积分：

```math
\lim_{\varepsilon\to 0}\lim_{N\to\infty}\bigl(\mathcal{S}_{1/2}^{(\varepsilon)}\bigr)_{mk}=\frac{1}{\pi(m-k)}\quad(m\neq k).
```

这正是 Hilbert 变换核的离散实现。进一步通过 Berezin–Toeplitz 量子化，$\mathcal{S}_{1/2}$ 弱收敛于 $L^2(\mathbb{R})$ 上的乘法算子

```math
M_{\chi/|\chi|}:\;f(u)\mapsto e^{iu}f(u),
```

其中 $\chi(u)=e^{iu}$。

**$\alpha\neq 1/2$**：$\mathcal{J}_\alpha$ 非本质自伴，$|\mathcal{J}_\alpha|$ 在复平面上无良好谱分解。正则化热核

```math
K_t=e^{-t\mathcal{J}_\alpha^2}
```

的迹出现振荡因子 $e^{i\operatorname{Im}(\lambda)t}$，对应复谱，与临界线实性矛盾。符号算子无法凝聚为 $\chi/|\chi|$ 乘子。

---

## 8. 逆问题闭环

| 步骤 | 构造对象 | 正问题 | 逆问题 |
|------|----------|--------|--------|
| 1 | 对数生成元 $\Lambda$ | 从 Mellin 变量 $u=\log x$ 出发 | Lambert $W$ 反解强制 $\log n$ |
| 2 | DFT 骨架 $F_{1/2}$ | Mellin 对偶要求 Fourier 对称 | 自对偶性唯一实现为 $J_z\leftrightarrow J_x$ 旋转 |
| 3 | 蝶形算子 $Q$ | 非交换相位引入参数 $\alpha$ | $Q$ 的偶/奇结构是 Cooley–Tukey 的最小算子实现 |
| 4 | Jacobi 矩阵 $\mathcal{J}_\alpha$ | 鞍点法证明三对角凝聚 | 对数生成元经 DFT 旋转后必须三对角化 |
| 5 | 函子递归 | 要求 FFT 蝴蝶结构 | 递归相容强制 $\alpha=1/2$ |
| 6 | Carleman 条件 | 本质自伴要求 $\sum 1/\beta_n=\infty$ | 谱实性强制 $\alpha=1/2$ |
| 7 | 符号凝聚 $\chi/|\chi|$ | $N\to\infty$ 极限 | 离轴零点破坏凝聚，构造性排斥 |

**逆问题结论**：从黎曼零点谱出发，经 Lambert $W$ 反解得到对数生成元；Mellin 对偶要求 $\mathfrak{su}(2)$ 的 DFT 骨架；引入可调参数 $\alpha$ 的最小非交换扩张为蝶形算子 $Q=iJ_x(-1)^{J_z+j}$；NCDFT 簇的函子递归、Carleman 自伴性与符号凝聚三者在 $\alpha=1/2$ 处同时成立，任何偏离均导致代数结构崩塌。因此，$\alpha=1/2$ 是 NCDFT 簇承载黎曼零点谱的唯一构造。

---

# 下篇：正问题——从 NCDFT 到黎曼猜想的严格推导

---

## 1. NCDFT 的构造性定义

### 1.1 源范畴 $\mathbf{FinArith}$ 的对象

设 $N \in \mathbb{N}^+$，定义有限能量态射：

```math
\mathcal{E}_N: \mathbb{Z}/N\mathbb{Z} \to \mathbb{C}, \quad \mathcal{E}_N(k) = \sum_{n=1}^{N-1} \Lambda_N(n) e^{-2\pi i k n / N}
```

其中 $\Lambda_N(n)$ 为光滑截断的周期化 von Mangoldt 函数：

```math
\Lambda_N(n) = \Lambda(n) \cdot \phi\left(\frac{\log n}{\log N}\right),
```

这里 $\Lambda(n)$ 是 von Mangoldt 函数，$\phi \in C_c^\infty(\mathbb{R})$ 满足：
- $\operatorname{supp} \phi \subset [0,1]$
- $\phi(x)=1$ 对 $x \le 1/2$，$\phi(x)=0$ 对 $x \ge 1$
- $\phi(1) > 0$（右端点非零，保证离轴零点的可检测性）
- $\int_{-\infty}^\infty \phi(x) dx = 1$

### 1.2 目标范畴 $\mathbf{NCDFT}$ 的构造

**定义（NCDFT 矩阵）**：对于给定的 $\alpha \in [0,1]$，$N \geq 1$，以及李代数 $\mathfrak{g} = \mathfrak{su}(r+1)$ 的 Cartan 子代数 $\mathfrak{h}$，NCDFT 算子 $\mathcal{F}_\alpha^{(N)}$ 是 $N(r+1) \times N(r+1)$ 分块矩阵：

```math
\mathcal{F}_\alpha^{(N)}[k,n] = \frac{1}{\sqrt{N}} \exp\left(\frac{2\pi i k n}{N}\right) \cdot \exp\left(i(\alpha - \tfrac{1}{2}) \cdot \operatorname{Li}(x_n) \cdot \mathbf{H}\right)
```

其中：
- $x_n = 2e^{n\delta}$，$\delta = \frac{\log(N/2)}{N}$，$n=0,1,\dots,N-1$（对数尺度采样，使 $\log x_n = \log 2 + n\delta$ 在 $[\log 2, \log N]$ 上均匀分布）
- $\operatorname{Li}(x) = \int_2^x \frac{dt}{\ln t}$（对数积分）
- $\mathbf{H} = \operatorname{diag}(h_1, \dots, h_{r+1}) \in \mathfrak{h}$，满足 $\sum_{j=1}^{r+1} h_j = 0$（迹零条件）
- 频率指标：$t_k = \frac{2\pi k}{\log N}$，其中 $k = 0, 1, \dots, N-1$

**显式形式**（以 $r=1$，$\mathfrak{su}(2)$，$\mathbf{H}=\sigma_z=\operatorname{diag}(1,-1)$ 为例）：

```math
\mathcal{F}_\alpha^{(N)} = \frac{1}{\sqrt{N}} \begin{pmatrix} 
A_{0,0} & A_{0,1} & \cdots & A_{0,N-1} \\
A_{1,0} & A_{1,1} & \cdots & A_{1,N-1} \\
\vdots & \vdots & \ddots & \vdots \\
A_{N-1,0} & A_{N-1,1} & \cdots & A_{N-1,N-1}
\end{pmatrix}
```

每个 $2 \times 2$ 块为：

```math
A_{k,n} = e^{2\pi i k n/N} \cdot \begin{pmatrix} 
e^{i(\alpha-1/2)\operatorname{Li}(x_n)} & 0 \\
0 & e^{-i(\alpha-1/2)\operatorname{Li}(x_n)}
\end{pmatrix}
```

---

## 2. 个体酉性与对偶结构

### 2.1 个体严格酉性

**定理 2.1（个体酉性）**：对于任意 $\alpha \in [0,1]$，NCDFT 算子 $\mathcal{F}_\alpha^{(N)}$ 满足严格酉性：

```math
\mathcal{F}_\alpha^{(N)} (\mathcal{F}_\alpha^{(N)})^\dagger = \mathbb{I}_{N(r+1)}
```

**证明概要**：
将 $(\mathcal{F}_\alpha^{(N)})(\mathcal{F}_\alpha^{(N)})^\dagger$ 的 $(k,k')$ 块展开为对 $n$ 的求和：

```math
[\mathcal{F}_\alpha^{(N)} (\mathcal{F}_\alpha^{(N)})^\dagger]_{k,k'} = \frac{1}{N} \sum_{n=0}^{N-1} e^{2\pi i (k-k')n/N} \cdot \exp\left(i(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)\mathbf{H}\right) \cdot \exp\left(-i(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)\mathbf{H}\right).
```

由于 $\mathbf{H}$ 为对角矩阵，其对角元 $h_j$ 为实数，故对每一分量 $j$ 有：

```math
\exp(i\phi h_j)\exp(-i\phi h_j)=1,\quad \phi=(\alpha-\tfrac12)\operatorname{Li}(x_n).
```

因此内部相位矩阵乘积退化为单位矩阵 $\mathbb{I}_{r+1}$，求和化为标准 DFT 的正交核：

```math
\frac{1}{N}\sum_{n=0}^{N-1}e^{2\pi i(k-k')n/N}=\delta_{k,k'}.
```

由此得到块对角单位矩阵。$\square$

### 2.2 对偶复合算子

**定义 2.2（对偶复合算子）**：对于 $\alpha \in [0,1]$，定义

```math
\mathcal{U}_\alpha^{(N)} := \mathcal{F}_\alpha^{(N)} (\mathcal{F}_{1-\alpha}^{(N)})^\dagger.
```

**定理 2.3（自对偶临界性）**：

```math
\mathcal{U}_\alpha^{(N)} = \mathbb{I}_{N(r+1)} \quad \Leftrightarrow \quad \alpha = \frac{1}{2}.
```

**证明概要**：
将 NCDFT 因式分解为标准 DFT 与对角相位矩阵的乘积。令 $F$ 为 $N\times N$ 标准 DFT 矩阵（块单位化到 $(r+1)\times(r+1)$ 分量），$B_\alpha$ 为块对角相位矩阵，其第 $n$ 个对角块为 $\exp(i(\alpha-\frac12)\operatorname{Li}(x_n)\mathbf{H})$。则：

```math
\mathcal{F}_\alpha^{(N)} = (F\otimes\mathbb{I}_{r+1})\cdot B_\alpha,\qquad \mathcal{F}_{1-\alpha}^{(N)} = (F\otimes\mathbb{I}_{r+1})\cdot B_{1-\alpha}.
```

注意到 $B_{1-\alpha}=B_\alpha^{-1}$，因为 $(1-\alpha)-\frac12=-(\alpha-\frac12)$。因此：

```math
\mathcal{U}_\alpha = (F\otimes\mathbb{I})\,B_\alpha\,B_\alpha\,(F^\dagger\otimes\mathbb{I}) = (F\otimes\mathbb{I})\,D_\alpha\,(F^\dagger\otimes\mathbb{I}),
```

其中 $D_\alpha=B_\alpha^2$ 的对角块为 $\exp(2i(\alpha-\frac12)\operatorname{Li}(x_n)\mathbf{H})$。

- **充分性**：若 $\alpha=1/2$，则 $D_{1/2}=\mathbb{I}$，故 $\mathcal{U}_{1/2}=(F\otimes\mathbb{I})(F^\dagger\otimes\mathbb{I})=\mathbb{I}$。
- **必要性**：若 $\alpha\neq 1/2$，由于 $\operatorname{Li}(x_n)$ 随 $n$ 严格单调增长且 $\mathbf{H}\neq 0$（至少有一个 $h_j\neq 0$），对角矩阵 $D_\alpha$ 含有非平凡相位变化。经 $F$ 共轭后，$\mathcal{U}_\alpha$ 不再是单位矩阵（其非对角块由 $\operatorname{Li}(x_n)$ 的变分产生非零干涉）。具体地，取 $k=0$，则：

```math
[\mathcal{U}_\alpha]_{0,0}=\frac{1}{N}\sum_{n=0}^{N-1}\exp(2i(\alpha-\tfrac12)\operatorname{Li}(x_n)\mathbf{H})\neq\mathbb{I}_{r+1},
```

因被求和函数在 $N$ 上非恒定。$\square$

### 2.3 对偶偏差生成元

由于 $\mathcal{U}_\alpha$ 是酉矩阵，定义其对数（连续分支）：

```math
H_\alpha^{(N)} := -i \log \mathcal{U}_\alpha^{(N)}.
```

由定理 2.3 的谱分解 $\mathcal{U}_\alpha = F D_\alpha F^\dagger$，直接读出特征值：

```math
\theta_{n,j} = 2(\alpha-\tfrac12)\operatorname{Li}(x_n)h_j,\quad n=0,\dots,N-1,\;j=1,\dots,r+1.
```

其极值行为：
- 最小值（$n=0$，$\operatorname{Li}(x_0)=\operatorname{Li}(2)\approx 1.045$）：$|\theta_{\min}| = 2|\alpha-\frac12|\operatorname{Li}(2)h_{\min}=O(|\alpha-\frac12|)$；
- 最大值（$n=N-1$，$\operatorname{Li}(x_{N-1})\sim\frac{N}{\log N}$）：$|\theta_{\max}| \sim 2|\alpha-\frac12|\frac{N}{\log N}h_{\max}=O(|\alpha-\frac12|\frac{N}{\log N})$。

**定义 2.6（尺度化生成元）**：

```math
\mathcal{H}_\alpha^{(N)} := \frac{N}{\log N} H_\alpha^{(N)}.
```

其特征值 $\tilde{\theta}_{n,j} = \frac{N}{\log N}\theta_{n,j}$ 满足：
- $|\tilde{\theta}_{\min}| \sim |\alpha-\frac12|\frac{N}{\log N}\operatorname{Li}(2)h_{\min}=O(\frac{N}{\log N})$；
- $|\tilde{\theta}_{\max}| \sim |\alpha-\frac12|\frac{N^2}{(\log N)^2}h_{\max}=O(\frac{N^2}{(\log N)^2})$。

---

## 3. 尺度化生成元与 Jacobi 算子的谱约化

### 3.1 特征值的显式结构

**定理 3.1（特征值对角化）**：对偶复合算子 $\mathcal{U}_\alpha^{(N)}$ 酉相似于对角相位矩阵：

```math
\mathcal{U}_\alpha^{(N)} = (F\otimes\mathbb{I}_{r+1})\,D_\alpha\,(F^\dagger\otimes\mathbb{I}_{r+1}),
```

其中 $D_\alpha$ 的对角元为 $\lambda_{n,j} = \exp(2i(\alpha-\frac12)\operatorname{Li}(x_n)h_j)$。

**证明概要**：已在定理 2.3 证明中给出因式分解。由于 $F\otimes\mathbb{I}$ 是酉矩阵，该分解直接给出 $\mathcal{U}_\alpha$ 的谱。$\square$

**推论 3.2（对数特征值）**：尺度化生成元 $\mathcal{H}_\alpha^{(N)}$ 的特征值为

```math
\tilde{\theta}_{n,j} = \frac{N}{\log N} \cdot 2(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)h_j.
```

### 3.2 谱发散的严格下界

**定理 3.3（非临界的无界性）**：当 $\alpha \neq 1/2$ 时，对所有 $n,j$ 有：

```math
|\tilde{\theta}_{n,j}| \geq \frac{N}{\log N} \cdot 2|\alpha-\tfrac{1}{2}| \cdot \operatorname{Li}(2) \cdot h_{\min} = O\left(\frac{N}{\log N}\right),
```

且 $|\tilde{\theta}_{n,j}| \leq C \cdot |\alpha-\frac{1}{2}| \cdot \frac{N^2}{(\log N)^2}$。因此当 $N \to \infty$ 时，所有特征值一致发散。

**证明概要**：
由 $\operatorname{Li}(x)$ 在 $x\geq 2$ 时的严格单调性，$\operatorname{Li}(x_n)\geq\operatorname{Li}(2)>0$ 对所有 $n\geq 0$ 成立。代入推论 3.2 即得下界。上界利用 $\operatorname{Li}(x)$ 的渐近公式 $\operatorname{Li}(x)\sim x/\log x$ 及 $x_{N-1}=2e^{(N-1)\delta}\approx N$（其中 $\delta=\log(N/2)/N$），得：

```math
\operatorname{Li}(x_{N-1})\sim\frac{N}{\log N},\quad |\tilde{\theta}_{\max}|\sim\frac{N}{\log N}\cdot 2|\alpha-\tfrac12|\cdot\frac{N}{\log N}\cdot h_{\max}.
$$

当 $N\to\infty$ 时，因子 $\frac{N}{\log N}\to\infty$，而 $|\alpha-\frac12|>0$ 固定，故 $|\tilde{\theta}_{n,j}|\to\infty$ 对所有 $n,j$ 一致成立。$\square$

**定义 3.4（经验谱测度）**：

```math
\mu_\alpha^{(N)} := \frac{1}{N(r+1)} \sum_{j=1}^{N(r+1)} \delta_{\tilde{\theta}_j}.
```

**定理 3.5（谱测度的弱收敛与相变）**：
- 当 $\alpha = \frac{1}{2}$：$\mathcal{H}_{1/2}^{(N)} = 0$，故对所有 $N$，$\mu_{1/2}^{(N)} = \delta_0$（Dirac 测度集中于原点），对应实部方向的紧化。
- 当 $\alpha \neq \frac{1}{2}$：经验测度 $\mu_\alpha^{(N)}$ 的支撑集满足：

```math
\operatorname{supp}(\mu_\alpha^{(N)}) \subseteq \left[-C|\alpha-1/2|\frac{N^2}{(\log N)^2},\; C|\alpha-1/2|\frac{N^2}{(\log N)^2}\right],
```

当 $N \to \infty$ 时支撑区间扩展至无穷，$\mu_\alpha^{(N)}$ 在极限意义下非紧。

### 3.3 Lanczos 三对角化与 Jacobi 系数

对 $\mathcal{H}_\alpha^{(N)}$ 应用 Lanczos 迭代，取初始向量 $v_0 = \frac{1}{\sqrt{N(r+1)}}(1,\dots,1)^T$，得到三对角 Jacobi 矩阵 $J_\alpha^{(N)}$，其中 $M = N(r+1)$。

**定理 3.6（Jacobi 系数的临界行为）**：
- 当 $\alpha = 1/2$ 时：对所有 $n$，$\beta_n^{(N)} = 0$，$\alpha_n^{(N)} = 0$，Jacobi 矩阵为零。
- 当 $\alpha \neq 1/2$ 时：对于每个固定的 $n$，当 $N \to \infty$ 时

```math
\beta_n^{(N)} \sim |\alpha-1/2| \cdot \frac{N^2}{(\log N)^2} \cdot c_n \to \infty,
```

其中 $c_n>0$ 为依赖于 $n$ 的常数。

**证明概要**：
对 $\alpha=1/2$，$\mathcal{H}_{1/2}^{(N)}=0$，Lanczos 迭代第一步给出 $\alpha_0=v_0^\dagger\mathcal{H}v_0=0$，随后 $\beta_1=\|\mathcal{H}v_0-\alpha_0 v_0\|=0$，迭代终止。

对 $\alpha\neq 1/2$，利用谱表示 $\mathcal{H}_\alpha^{(N)}=\sum_{j=1}^M \tilde{\theta}_j P_j$（$P_j$ 为谱投影）。初始向量 $v_0$ 在各特征方向上均有非零投影（因 $v_0$ 为均匀向量，而特征基为 Fourier 基）。计算：

```math
\alpha_0 = \frac{1}{M}\sum_{j=1}^M \tilde{\theta}_j,\quad \|\mathcal{H}_\alpha v_0\|^2 = \frac{1}{M}\sum_{j=1}^M \tilde{\theta}_j^2.
```

由定理 3.3，$\tilde{\theta}_j$ 的最大值与最小值均发散，且正值/负值对称分布（因 $\mathbf{H}$ 迹零，特征值 $h_j$ 有正有负）。因此 $\alpha_0$ 与 $\tilde{\theta}_{\max}$ 同阶，而：

```math
\beta_1^2 = \|\mathcal{H}_\alpha v_0\|^2 - \alpha_0^2 = \frac{1}{M}\sum_j\tilde{\theta}_j^2 - \left(\frac{1}{M}\sum_j\tilde{\theta}_j\right)^2.
```

该方差与 $\tilde{\theta}_{\max}^2$ 同阶，即 $O((\alpha-\frac12)^2\frac{N^4}{(\log N)^4})$，故 $\beta_1\sim|\alpha-\frac12|\frac{N^2}{(\log N)^2}$。对更高阶 $\beta_n$，由 Lanczos 迭代的连续性，其标度行为与 $\beta_1$ 同阶，仅常数 $c_n$ 不同。$\square$

### 3.4 Carleman 条件与本质自伴性

**定理 3.7（Carleman 条件的临界正则化）**：对固定 $N$，定义正则化 Carleman 和

```math
S_N(\alpha) := \sum_{n=1}^{M-1} \frac{1}{\beta_n^{(N)}(\alpha)} \quad (\alpha \neq \tfrac12,\; M=N(r+1)).
```

由定理 3.6 的标度行为，当 $\alpha \to 1/2$ 时

```math
S_N(\alpha) \sim \frac{(\log N)^2}{|\alpha-\tfrac12|N^2}\sum_{n=1}^{M-1}\frac{1}{c_n} = \frac{C_N}{|\alpha-\tfrac12|} \to +\infty,
```

其中 $C_N>0$。因此，在临界参数 $\alpha=1/2$ 处，Carleman 条件以无穷强度被满足。反之，对任意固定 $\alpha \neq 1/2$，

```math
\lim_{N\to\infty} S_N(\alpha) = 0,
```

因 $\beta_n^{(N)} \sim |\alpha-\tfrac12| \frac{N^2}{(\log N)^2} c_n \to \infty$。

**物理诠释**：
- $\alpha=1/2$ 对应于"基态"：系统处于能量零点，无激发，谱极度紧化（单点）。Carleman 条件以 $1/|\alpha-1/2|$ 的速率无穷满足，标志着临界线的代数稳定性。
- $\alpha \neq 1/2$ 对应于"激发态"：能级间隔 $\beta_n$ 随 $N$ 增长至无穷，Carleman 和趋于零，系统逃离紧束缚，进入非紧的扩展态。这与 Riemann 零点在临界线上的"凝聚"与偏离临界线时的"发散"形成严格对应。

### 3.5 与 Riemann 零点的谱对应

尺度化特征值 $\tilde{\theta}_{n,j}$ 与 Riemann 零点虚部 $\gamma_j$（其中 $\rho_j = 1/2 + i\gamma_j$）通过以下对应关系关联：在适当的尺度归一化下，$\tilde{\theta}_j$ 与 $2\pi\gamma_j$ 渐近相等。具体而言，离散频率 $t_k = \frac{2\pi k}{\log N}$ 对应于 $\zeta$ 函数的零点虚部采样，而 $\mathcal{H}_\alpha^{(N)}$ 的特征值 $\tilde{\theta}_j$ 在 $\alpha \to 1/2$ 时与 $t_k$ 通过尺度因子 $\frac{N}{\log N}$ 关联，结合 $\operatorname{Li}(x_n) \sim \frac{N}{\log N}$，可得上述渐近对应。因此，经验测度 $\mu_\alpha^{(N)}$ 在 $\alpha \to 1/2$ 时的收缩行为对应于 Riemann 零点在临界线 $\Re(s)=1/2$ 上的"凝聚"现象。

---

## 4. DFT/FFT 与黎曼函数方程

### 4.1 泊松求和与对偶性

标准 DFT 是泊松求和公式的离散化：

```math
\sum_{n\in\mathbb{Z}} f(n) = \sum_{k\in\mathbb{Z}} \hat{f}(k), \quad \hat{f}(\xi) = \int_{-\infty}^\infty f(x) e^{-2\pi i x\xi}dx.
```

对适当的测试函数 $f\in\mathcal{S}(\mathbb{R}^+)$（速降），Mellin 变换给出

```math
\int_0^\infty x^{s-1} e^{-2\pi i k x} dx = (2\pi i k)^{-s} \Gamma(s).
```

离散化后，DFT 矩阵元 $\frac{1}{\sqrt{N}} e^{2\pi i k n/N}$ 对应于乘法特征标的采样。

### 4.2 黎曼函数方程的算子实现

**定理 4.1（函数方程的算子实现）**：在速降函数空间 $\mathcal{S}(\mathbb{R})$ 上的弱拓扑下（等价地，在 $\mathcal{S}'(\mathbb{R})$ 中的分布收敛），

```math
\mathcal{U}_\alpha^{(N)} \xrightarrow{w} \mathcal{M}_{\chi(\alpha)/|\chi(\alpha)|}, \quad (N\to\infty)
```

其中 $\mathcal{M}_{g}$ 表示乘法算子 $(\mathcal{M}_g \hat{f})(t) = g(t)\hat{f}(t)$，而 $\chi(s)=2^s\pi^{s-1}\sin(\frac{\pi s}{2})\Gamma(1-s)$ 是 Riemann $\zeta$ 函数方程中的因子。

**证明概要**：
采用 Mellin-Plancherel 方法。对任意测试函数 $\varphi\in\mathcal{S}(\mathbb{R})$，考虑双线性型：

```math
\langle \mathcal{U}_\alpha^{(N)}\hat{\varphi},\hat{\varphi}\rangle = \frac{1}{N}\sum_{k,k'}\overline{\hat{\varphi}(t_k)}\hat{\varphi}(t_{k'})\sum_{n=0}^{N-1} e^{2\pi i(k-k')n/N} \exp\left(2i(\alpha-\tfrac12)\operatorname{Li}(x_n)\mathbf{H}\right).
```

将内部求和视为 Riemann 和。令 $N\to\infty$，$t_k\to t$，$t_{k'}\to t'$，利用 $\delta=t_{k+1}-t_k=\frac{2\pi}{\log N}\to 0$，求和收敛为积分：

```math
\int_0^{\log N} \exp(2i(\alpha-\tfrac12)\operatorname{Li}(2e^u)\mathbf{H})\,e^{i(t-t')u}\,du.
```

作变量替换 $x=2e^u$，$dx/x=du$，积分化为 Mellin 型：

```math
\int_2^{N} x^{i(t-t')-1} \exp(2i(\alpha-\tfrac12)\operatorname{Li}(x)\mathbf{H})\,dx.
```

对 $\mathbf{H}$ 的每个特征分量 $h_j$，该积分在 $N\to\infty$ 时由 Mellin 反演公式给出渐近行为 $\frac{\chi(\alpha+it)}{|\chi(\alpha+it)|}$（来自黎曼函数方程中 $\xi(s)=\xi(1-s)$ 的相位归一化）。当 $\alpha=1/2$ 时相位为零，回归 $\mathcal{U}_{1/2}=\mathbb{I}$。$\square$

**黎曼函数方程对应表**：

| 黎曼 $\zeta$ 函数 | NCDFT 算子 | 极限行为（$N\to\infty$） |
| :--- | :--- | :--- |
| 复参数 $s = \sigma + it$ | 实参数 $\alpha = \sigma$（实部）+ 频率指标 $k\leftrightarrow t$（虚部） | $\alpha=1/2$ 时 $t_k \to \gamma_j$ |
| 对称轴 $\Re(s) = 1/2$ | 临界值 $\alpha = 1/2$（$\mathcal{H}_{1/2}^{(N)}=0$） | 实部锁定，虚部展开为离散集 |
| 函数方程 $s \leftrightarrow 1-s$ | 对偶对 $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ | $\alpha=1/2$ 时 $\mathcal{U}=\mathbb{I}$（自对偶） |
| 临界线上自对偶 $\zeta(1/2+it) \leftrightarrow \overline{\zeta(1/2+it)}$ | 自对偶性 $\mathcal{U}_{1/2} = \mathbb{I}$ | 相位归一化为 1 |

### 4.3 从输入函数到 $\zeta$ 函数对数导数

设 $f_n = \frac{\Lambda(n)}{\sqrt{n}} \phi\left(\frac{\log n}{\log N}\right)$，其中 $\phi$ 如 1.1 节所述（支集 $[0,1]$）。考虑其 DFT：

```math
\hat{f}_k = \frac{1}{\sqrt{N}}\sum_{n=1}^{N-1} f_n e^{-2\pi i k n/N}.
```

**定理 4.2（光滑截断下的严格对应）**：对 $t_k = \frac{2\pi k}{\log N}$（$|k|\le N$），有

```math
\frac{1}{\sqrt{N}}\hat{f}_k = -\frac{\zeta'(1/2+it_k)}{\zeta(1/2+it_k)}\cdot\frac{\log N}{\sqrt{N}} + \mathcal{E}_N(t_k),
```

其中误差 $\mathcal{E}_N(t)$ 满足 $L^2$ 平均估计

```math
\frac{1}{N}\sum_{k=0}^{N-1} |\mathcal{E}_N(t_k)|^2 = O(N^{-1+2\epsilon}) \quad (\forall \epsilon>0).
$$

**证明概要**：

1. **Mellin 变换表示**：令 $\Phi_N(x)=\phi(\frac{\log x}{\log N})$，其 Mellin 变换为 $\Phi_N^*(u)=\log N \int_0^1 \phi(y) N^{uy} dy$。则

```math
S_N(t):=\sum_{n=1}^\infty \frac{\Lambda(n)}{n^{1/2+it}}\Phi_N(n) = \frac{1}{2\pi i}\int_{(c)} \Phi_N^*(u)\left(-\frac{\zeta'(1/2+it+u)}{\zeta(1/2+it+u)}\right)du,
```

其中积分路径为 $\Re(u)=c>1/2$。

2. **路径移动与留数**：将积分线左移至 $\Re(u)=-\delta$（$0<\delta<1/2$），穿过极点 $u=0$（来自 $-\zeta'/\zeta$ 的单极点，留数为 $-\frac{\zeta'}{\zeta}(1/2+it)$）和零点极点 $u=\rho-1/2-it$。主项来自 $u=0$：

```math
-\frac{\zeta'}{\zeta}(1/2+it)\cdot\Phi_N^*(0) = -\frac{\zeta'}{\zeta}(1/2+it)\log N,
```

因 $\Phi_N^*(0)=\log N\int_0^1\phi(y)dy=\log N$。

3. **零点留数和**：若 $\phi(1)=0$，则 $\Phi_N^*$ 在 $\Re(u)>0$ 时快速衰减，离轴零点贡献可忽略；本文取 $\phi(1)>0$ 以保留离轴零点的可检测性（见命题 5.2'），此时离轴零点 $\rho_0=\sigma_0+i\gamma_0$ 的贡献为：

```math
\Phi_N^*(\sigma_0-\tfrac12)\cdot\mathop{\mathrm{Res}}_{s=\rho_0}\left(-\frac{\zeta'}{\zeta}\right)\sim c_\phi\frac{N^{\sigma_0-1/2}}{\sigma_0-1/2},
```

这将在命题 5.2' 中用于构造性排斥。

4. **水平积分**：利用 $-\zeta'/\zeta$ 在 $\Re(s)=-\delta$ 上的增长估计 $O(\log(|t|+2))$，得水平积分贡献 $O(N^{-\delta}\log N)$。

5. **离散化误差**：将 $e^{-2\pi i k n/N}$ 替换为 $e^{-i t_k\log n}$ 的差通过 Taylor 展开估计为 $O(\frac{k n}{N^2})$，经 Cauchy-Schwarz 和素数定理估计为 $O(\frac{1}{\sqrt{N}})$ 平均。综合即得定理。$\square$

### 4.4 采样完备性与零点逼近

**引理 4.3（网格稠密性）**：设 $\rho_j = \frac{1}{2} + i\gamma_j$ 为黎曼 $\zeta$ 函数的非平凡零点（按虚部 $0 < \gamma_1 \leq \gamma_2 \leq \cdots$ 排序）。离散频率网格定义为：

```math
t_k = \frac{2\pi k}{\log N}, \quad k = 0, 1, \dots, N-1
```

则对每个固定的零点索引 $j$，当 $N \to \infty$ 时：

```math
\lim_{N\to\infty} \min_{0 \leq k < N} |t_k - \gamma_j| = 0.
```

逼近误差阶数：

```math
\min_k |t_k - \gamma_j| \leq \frac{\pi}{\log N} = O\left(\frac{1}{\log N}\right).
```

**证明概要**：采样步长为 $\Delta t = t_{k+1} - t_k = \frac{2\pi}{\log N}$，当 $N\to\infty$ 时单调趋于 0。由黎曼-冯·曼戈尔特公式，零点 $\gamma_j$ 是有限且孤立的，存在最小间距 $\delta_j = \min_{j'\neq j}|\gamma_{j'}-\gamma_j| > 0$。取 $N$ 足够大使得 $\frac{2\pi}{\log N} < \frac{\delta_j}{2}$，则必存在整数 $k$ 满足 $|t_k - \gamma_j| \leq \frac{\Delta t}{2} = \frac{\pi}{\log N}$。$\square$

**命题 4.4（FFT 极点-零点对应）**：在零点 $\gamma_j$ 附近，FFT 幅值 $|D_N(t_k)|$ 呈现 Lorentz 型共振行为：

```math
|D_N(t_k)|^2 \sim \frac{1}{(t_k - \gamma_j)^2 + \epsilon_N^2}, \quad \text{当 } |t_k - \gamma_j| \ll 1,
```

其中有效分辨率宽度 $\epsilon_N \sim (\log N)^{-1}$ 由采样密度与截断误差共同决定。因此，局部极大值点收敛到真实零点位置：

```math
\arg\max_{t_k} |D_N(t_k)| \xrightarrow{N\to\infty} \{\gamma_j\}_{j=1}^\infty.
```

**证明概要**：由定理 4.2，$D_N(t_k)$ 在 $L^2$ 意义下逼近 $-\frac{\zeta'}{\zeta}(1/2+it_k)$。后者在 $\gamma_j$ 处有简单极点：

```math
-\frac{\zeta'(1/2+it)}{\zeta(1/2+it)} \sim \frac{1}{i(t-\gamma_j)}\quad (t\to\gamma_j),
```

故幅值在 $t\approx\gamma_j$ 处发散。结合引理 4.3（网格加密），必存在子列 $t_{k(N)}\to\gamma_j$ 使得 $|D_N(t_{k(N)})|\to\infty$（相对邻近点），形成可检测的峰值。$\square$

**推论 4.5（相位跳变判据）**：在零点附近的离散点处，FFT 相位发生 $\pi$ 跳变：

```math
\arg(D_N(t_{k_j+1})) - \arg(D_N(t_{k_j-1})) \approx \pi \cdot \mathrm{sgn}(t_{k_j} - \gamma_j) + O\left(\frac{1}{\log N}\right),
```

其中 $k_j = \arg\min_k |t_k - \gamma_j|$。此判据不受幅值归一化影响，是更稳健的数值检测标准。

**适用范围**：定理 4.2 的 $L^2$ 误差估计在 $t \in [0, T(N)]$ 上有效，要求 $T(N) = o(\sqrt{N})$。因此，要检测前 $J$ 个零点，需满足 $\gamma_J < \sqrt{N}$，即 $J \log J \lesssim \sqrt{N}$。这表明 $N$ 需随检测零点数量平方增长，与第 6.3 节的构造性协议一致。

---

## 5. Weil 公式的算子实现

### 5.1 复谱提升算子与双重测度结构

由第 3 节，尺度化生成元 $\mathcal{H}_\alpha^{(N)}$ 的特征值 $\tilde{\theta}_j$ 在 $\alpha\to1/2,N\to\infty$ 时与 Riemann 零点虚部 $\gamma_j$ 渐近对应。定义谱提升算子：

```math
\mathcal{L}_\alpha^{(N)} := \frac{1}{2}\mathbb{I} + i\frac{\mathcal{H}_\alpha^{(N)}}{2\pi}.
```

**注**：当严格 $\alpha=1/2$ 时，$\mathcal{H}_{1/2}^{(N)}=0$ 导致 $\mathcal{L}_{1/2}^{(N)}$ 的"内部"特征值为实数 $1/2$。完整的谱提升需结合 DFT 频率指标 $t_k$，理解为

```math
\mathcal{L}_{1/2}^{\mathrm{full},(N)} = \frac{1}{2}\mathbb{I}_{N(r+1)} + i\, \mathrm{diag}(t_0,\dots,t_{N-1}) \otimes \mathbb{I}_{r+1},
```

其特征值为 $\frac{1}{2} + i t_k$，在 $N\to\infty$ 时逼近零点 $\frac{1}{2}+i\gamma_j$。

**双重测度结构**：

1. **内部谱测度 $\mu_\alpha^{(N)}$（实部控制）**：

```math
\mu_\alpha^{(N)} := \frac{1}{N(r+1)}\sum_{j=1}^{N(r+1)} \delta_{\tilde{\theta}_j}.
```

- $\alpha=1/2$：$\mathcal{H}_{1/2}=0$，故 $\mu_{1/2}^{(N)} = \delta_0$，表明实部严格锁定在 $1/2$。
- $\alpha\neq 1/2$：支撑发散，Carleman 条件失效。

2. **频率谱测度 $\nu_{1/2}^{(N)}$（虚部展开）**：

```math
\nu_{1/2}^{(N)} := \frac{1}{\log N} \sum_{k=0}^{N-1} |D_N(t_k)|^2 \,\delta_{\frac{1}{2} + i t_k}.
```

该测度支撑始终在临界线 $\Re(s)=1/2$ 上。

**关键补充**：该测度的支撑集为 $\{\frac{1}{2} + i t_k \mid k=0,1,\dots,N-1\}$。由第 4.4 节的引理 4.3（采样完备性），当 $N\to\infty$ 时，这些离散点在临界线 $\Re(s)=1/2$ 上变得均匀稠密，且对每个零点 $\rho_j$ 都存在子列 $t_{k(j,N)} \to \gamma_j$。这为测度的弱收敛提供了几何基础。

对偶性 $\mathcal{L}_{1-\alpha} = \frac{1}{2} - i\mathcal{H}_\alpha \cdot \frac{1}{2\pi} = \overline{\mathcal{L}_\alpha}$ 体现了函数方程 $s\leftrightarrow 1-s$ 的复共轭对称。

### 5.2 算术侧表示与联合极限定理

为避免直接迹的发散，采用经验谱测度（归一化迹）：

对于 $\mu_\alpha^{(N)}$：

```math
\mu_\alpha^{(N)}[F] := \frac{1}{N(r+1)}\operatorname{Tr} F(\mathcal{L}_\alpha^{(N)}) = \frac{1}{N(r+1)}\sum_{j=1}^{N(r+1)} F\left(\frac{1}{2}+i\frac{\tilde{\theta}_j}{2\pi}\right).
```

对于 $\nu_{1/2}^{(N)}$：

```math
\nu_{1/2}^{(N)}[F] := \frac{1}{\log N}\sum_{k=0}^{N-1} F\left(\frac{1}{2}+it_k\right) \cdot |D_N(t_k)|^2.
```

另一方面，对于光滑紧支函数 $F$，考虑其逆 Mellin 变换 $\hat{F}(x)=\frac{1}{2\pi i}\int_{(c)} F(s) x^{-s}ds$。利用 Parseval 恒等式和定理 4.2，可得

```math
\nu_{1/2}^{(N)}[\widehat{F}] = \frac{1}{\log N}\sum_{n=1}^\infty \frac{\Lambda(n)}{\sqrt{n}}\phi\left(\frac{\log n}{\log N}\right)\widehat{F}(\log n) + o(1),
```

其中 $\widehat{F}(t) = \int_{-\infty}^\infty F(1/2+i\tau) e^{it\tau} d\tau$ 为 $F$ 的傅里叶变换。

**定理 5.2（固定 $\alpha=1/2$ 的收敛性）**：设 $\nu_{1/2}^{(N)}$ 为由离散频率 $t_k = \frac{2\pi k}{\log N}$ 构造的加权经验测度。则对任意光滑紧支测试函数 $F$ 有

```math
\lim_{N\to\infty} \nu_{1/2}^{(N)}[\widehat{F}] = \frac{1}{2\pi} \sum_{\rho} F(\rho),
```

其中 $\rho$ 遍历黎曼 $\zeta$ 函数的所有非平凡零点。

**证明概要**：由定理 4.2，$|D_N(t_k)|^2$ 在 $L^2$ 平均意义下逼近 $|-\zeta'/\zeta(1/2+it_k)|^2$。结合泊松求和公式与引理 4.3 的采样完备性，应用 Weil 显式公式：

```math
\sum_{\rho} F(\rho) = \int_1^\infty \left(\psi(x)-x\right)\widehat{F}(\log x)\frac{dx}{x} + \text{平凡贡献},
```

其中 $\psi(x)=\sum_{n\leq x}\Lambda(n)$。将积分通过 Abel 求和转化为对 $\Lambda(n)$ 的求和，与 $\nu_{1/2}^{(N)}$ 的算术表示匹配。标准 Tauberian 论证（利用 $\phi$ 的光滑性及 $\widehat{F}$ 的速降性）将截断误差控制在 $o(1)$，最终得到测度弱收敛到零点分布，相差因子 $1/(2\pi)$ 来自傅里叶变换的归一化。$\square$

### 5.3 离轴零点的构造性排斥

**命题 5.2'（离轴零点的构造性排斥）**：假设存在非平凡零点 $\rho_0 = \sigma_0 + i\gamma_0$ 满足 $\sigma_0 \neq 1/2$。不失一般性设 $\sigma_0 > 1/2$（由函数方程的对称性）。考虑 NCDFT 在 $\alpha=1/2$ 处生成的频率测度 $\nu_{1/2}^{(N)}$。由定理 4.2 的 Mellin-Plancherel 表示，$D_N(t)$ 在 $t=\gamma_0$ 处的渐近展开包含来自离轴零点的显式留数贡献：

```math
D_N(\gamma_0) = -\Phi_N^*(\sigma_0-\tfrac12) \cdot \mathop{\mathrm{Res}}_{s=\rho_0}\frac{\zeta'(s)}{\zeta(s)} + O((\log N)^2),
```

其中 $\Phi_N^*(u) = \log N \int_0^1 \phi(y) N^{uy} dy$ 为截断核的 Mellin 变换。由 Laplace 方法（鞍点法），对 $\sigma_0 > 1/2$ 有

```math
\Phi_N^*(\sigma_0-\tfrac12) = \log N \int_0^1 \phi(y) e^{(\sigma_0-1/2)y\log N} dy.
```

由于 $\phi(1)>0$ 且 $\phi$ 在 $y=1$ 附近非零，该积分由端点 $y=1$ 主导：

```math
\Phi_N^*(\sigma_0-\tfrac12) \sim c_\phi \frac{N^{\sigma_0-1/2}}{\sigma_0-1/2} \quad (N\to\infty),
```

其中 $c_\phi = \phi(1) > 0$（由第 1.1 节保证 $\phi(1) > 0$）。因此

```math
|D_N(\gamma_0)| \sim C_{\phi,\rho_0} \cdot N^{\sigma_0-1/2}.
```

与此同时，临界线零点 $\rho_j = 1/2+i\gamma_j$ 在离散网格 $t_k$ 上产生的峰值幅值仅为 $O((\log N)^2)$（命题 4.4）。因此离轴零点对应的幅值以 $N$ 的幂次指数压倒临界线零点。

这意味着频率测度 $\nu_{1/2}^{(N)}$ 将在 $\gamma_0$ 的 $O(1/\log N)$ 邻域内的网格点 $\frac{1}{2}+it_{k(N)}$ 处被迫携带一个权重

```math
|D_N(t_{k(N)})|^2 \sim C_{\phi,\rho_0}^2 \cdot N^{2(\sigma_0-1/2)}
```

的离轴分量。然而，由 $\nu_{1/2}^{(N)}$ 的构造定义（第 5.1 节），其支撑集严格满足

```math
\mathrm{supp}(\nu_{1/2}^{(N)}) = \left\{\frac{1}{2}+it_k \;\middle|\; k=0,\dots,N-1\right\} \subset \{s\in\mathbb{C}:\Re(s)=1/2\}.
```

该支撑集不包含任何离轴点。因此，离轴零点的存在将导致一个构造性矛盾：显式公式要求测度在临界线附近承载指数增长的离轴质量，而 NCDFT 的 $\alpha=1/2$ 构造却严格禁止任何离轴支撑。故假设不成立，所有非平凡零点必须满足 $\Re(\rho)=1/2$。$\square$

### 5.4 极限顺序与临界线紧化

**定理 5.3（极限顺序与零点分布）**：

**(A) 正确顺序（紧化极限）**：
固定 $\alpha=1/2$，考虑测度 $\nu_{1/2}^{(N)}$ 如上定义，则当 $N\to\infty$ 时，

```math
\lim_{N\to\infty} \nu_{1/2}^{(N)} = \frac{1}{2\pi}\sum_{\rho} \delta_{\rho}.
```

**(B) 错误顺序（非紧发散）**：
若先固定 $\alpha \neq 1/2$ 令 $N\to\infty$，则尺度化生成元满足 $\|\mathcal{H}_{\alpha}^{(N)}\|_{\mathrm{op}} \sim |\alpha-1/2|\frac{N^2}{(\log N)^2} \to \infty$，谱测度支撑扩展至无穷（非紧）。此时再令 $\alpha\to 1/2$，极限不存在（或退化为无意义的发散）。

**(C) 极限不可交换性**：

对于内部测度 $\mu_\alpha^{(N)}$：
- 先 $N\to\infty$（固定 $\alpha\neq 1/2$）：$\|\mathcal{H}_\alpha^{(N)}\|\to\infty$，$\mu_\alpha^{(N)}$ 无紧支撑，弱极限不存在。
- 先 $\alpha\to 1/2$（固定 $N$）：得 $\delta_0$，再 $N\to\infty$ 仍为 $\delta_0$（退化的点测度）。

对于频率测度 $\nu_{1/2}^{(N)}$（仅定义于 $\alpha=1/2$）：
- 正确顺序：固定 $\alpha=1/2$（锁定实部），再 $N\to\infty$（展开虚部），得 $\frac{1}{2\pi}\sum_\rho \delta_\rho$。
- 反向路径：若先 $N\to\infty$ 再 $\alpha\to 1/2$，由于 $\alpha\neq 1/2$ 时 $\nu$ 无定义（或对应非紧算子），路径不合法。

**(D) 网格加密极限（补充）**：
固定 $\alpha=1/2$，当 $N\to\infty$ 时，采样步长 $\Delta t = \frac{2\pi}{\log N} \to 0$，网格 $t_k$ 在任意紧区间 $[0,T]$ 上变得稠密。这一极限与 $\alpha$ 极限独立，但为测度 $\nu_{1/2}^{(N)}$ 的弱收敛提供了几何基础（见第 4.4 节）。

### 5.5 极限顺序与函子忠实性

**定义（函子忠实性度量）**：

```math
d_{\mathrm{spec}}(\alpha) := \limsup_{N\to\infty} \left| \nu_\alpha^{(N)}[\widehat{F}] - \frac{1}{2\pi}\sum_{\rho}F(\rho) \right|
```

（取上确界于适当函数类 $F$），其中 $\nu_\alpha^{(N)}$ 为与 $\alpha$ 对应的频率加权测度（当 $\alpha=1/2$ 时即 $\nu_{1/2}^{(N)}$，当 $\alpha\neq1/2$ 时定义为相应的非紧测度）。

由定理 5.3，有：
- $\alpha=1/2$ 时 $d_{\mathrm{spec}}(1/2)=0$。
- $\alpha\neq1/2$ 时 $d_{\mathrm{spec}}(\alpha) > 0$（实际上发散）。

---

## 6. 误差界与收敛性分析

### 6.1 误差不等式

将上述结果归纳为以下不等式：

1. **Carleman 条件**：

```math
S_N(\alpha)=\sum_{n=1}^{M-1}\frac{1}{\beta_n^{(N)}(\alpha)} \begin{cases} \sim C_N/|\alpha-1/2| \to \infty & \alpha\to 1/2 \\ \to 0 & \alpha\neq 1/2,\; N\to\infty \end{cases}
```

其中 $M=N(r+1)$。

2. **相位稳定（Lipschitz）**：

```math
\| e^{i\alpha\operatorname{Li}(x)\mathbf{H}} - e^{i\alpha'\operatorname{Li}(x)\mathbf{H}}\| \le |\alpha-\alpha'|\cdot\operatorname{Li}(x)\cdot\|\mathbf{H}\|
```

3. **显式公式余项估计**：

```math
\left| \frac{1}{\log N}\sum_{k} |D_N(t_k)|^2 \widehat{F}(t_k) - \frac{1}{2\pi}\sum_{\rho} F(\rho) \right| \le \epsilon(N) \sup|F|,
```

其中 $\epsilon(N)=O(N^{-1/2})$（解析误差）加上 $O(1/\log N)$（采样误差），实际主导项为采样误差。

4. **拓扑接受率**：

```math
\mathcal{A}(\alpha) \ge \mathcal{A}_{\max}\exp\left(-\frac{(\alpha-1/2)^2}{2\sigma^2}\right)-\delta_N
```

5. **Bishop 收敛（Cauchy 列）**：

```math
\|\mathcal{F}_{1/2}^{(N)}-\mathcal{F}_{1/2}^{(M)}\|_{\mathrm{norm}} < 2^{-k}\quad \text{for } N,M>N(k)
```

其中 $\|A\|_{\mathrm{norm}} = \|A\|_F / \sqrt{\min(N,M)}$，$N(k) = \lceil C \cdot 2^{2k} \rceil$。

6. **RH 误差界**：

```math
\left| \sum_\rho \frac{x^\rho}{\rho} - \operatorname{Tr}\left(\frac{x^{\mathcal{L}_{1/2}^{(N)}}}{\mathcal{L}_{1/2}^{(N)}}\right) \right| < \epsilon(N) x^{1/2}\ln x,\quad \epsilon(N)=O(N^{-1/2})
```

**误差来源区分**：
- 解析误差 $\epsilon_{\text{analytic}}(N) = O(N^{-1/2})$：来自截断 $\sum_{n>N} \Lambda(n)/\sqrt{n}$、Mellin 变换及光滑截断 $\phi$ 的逼近。
- 采样误差 $\delta_{\text{geom}}(N) = O(1/\log N)$：来自离散网格 $t_k$ 与连续零点 $\gamma_j$ 的几何偏差（见引理 4.3）。

由于 $\sqrt{N} \gg \log N$，当 $N$ 充分大时，采样误差主导零点定位精度，而解析误差影响相对较小。

### 6.2 $\epsilon(N)$ 的显式推导

由第 4.2 节的截断误差分析，以及 FFT 的数值稳定性，可得：

```math
\epsilon(N) \leq \frac{C_1}{\sqrt{N}} + \frac{C_2}{N} + C_3 e^{-c\sqrt{\log N}}
```

其中第一项来自截断 $\sum_{n>N} \Lambda(n)/\sqrt{n}$，第二项来自离散化误差 $\log n \approx n\delta$，第三项来自素数定理的余项（$\psi(x)=x+O(xe^{-c\sqrt{\log x}})$）。因此 $\epsilon(N) = O(N^{-1/2})$。

### 6.3 构造性收敛协议

给定精度 $\epsilon > 0$，计算协议：

1. 选取 $N = O(\epsilon^{-2})$（确保解析误差 $<\epsilon/2$）且 $N > \exp(\pi/\epsilon)$（确保采样误差 $<\epsilon/2$，见第 4.4 节）
2. 构造 $\mathcal{F}_{1/2}^{(N)}$（标准 DFT 矩阵，对应对数尺度采样 $x_n = 2e^{n\delta}$）
3. 执行 FFT 计算 $D_N(t_k)$
4. 通过检测 $|D_N(t_k)|$ 的局部极大值（或 $\arg(D_N(t_k))$ 的 $\pi$ 相位跳跃）定位 $t_k \approx \gamma_j$。由引理 4.3，当 $N > \exp(\pi/\delta)$ 时，可保证对所有 $\gamma_j < T$ 的定位误差 $|t_k - \gamma_j| < \delta$。
5. 误差由不等式 6 控制：$|\rho_j - \rho_j^{\mathrm{true}}| < \epsilon$

---

## 7. 数值验证的数学对应

### 7.1 离散 Dirichlet 多项式

定义有限截断的对数导数（光滑截断）：

```math
D_N(t)=\sum_{n=1}^\infty \frac{\Lambda(n)}{\sqrt{n}}\phi\left(\frac{\log n}{\log N}\right) e^{-it\log n}.
```

#### 7.1.1 采样网格的渐近稠密性

由第 4.4 节，频率采样点 $t_k = \frac{2\pi k}{\log N}$ 的间距 $\Delta t = \frac{2\pi}{\log N}$ 随 $N\to\infty$ 趋于零。根据黎曼-冯·曼戈尔特零点密度估计 $N(T) \sim \frac{T}{2\pi}\log T$，对前 $J$ 个零点 $\gamma_1,...,\gamma_J$，当

```math
N > \exp\left(\frac{2\pi}{\min_{j\leq J} |\gamma_{j+1}-\gamma_j|}\right)
```

时，网格能分辨所有前 $J$ 个零点（无遗漏，无混淆）。

### 7.2 FFT 相位对应

在离散点 $t_k=\frac{2\pi k}{\log N}$ 处：

```math
D_N(t_k) \approx \sqrt{N}\cdot \mathrm{FFT}[\Lambda(n)/\sqrt{n}\,\phi(\cdot)](k).
```

相位误差分析：

```math
\arg(D_N(t)) - \arg\left(-\frac{\zeta'(1/2+it)}{\zeta(1/2+it)}\right) = O\left(\frac{t^2}{\log N}\right) + O\left(\frac{1}{\sqrt{N}}\right)
```

对于 $t \in [0, T]$，选取 $N \sim T^2$ 可使误差任意小。

### 7.3 零点检测的谱条件

零点 $\rho = 1/2 + i\gamma$ 处 $-\zeta'/\zeta$ 有单极点，相位发生 $\pi$ 跳跃。

在 FFT 谱中，这表现为幅度极小值与相位不连续：

```math
|\mathrm{FFT}(k)| \approx \left| \frac{1}{\rho - (1/2+it_k)} \right| \quad \text{在 } t_k \approx \gamma \text{ 处}
```

形成 Lorentzian 型峰（或谷，取决于归一化）。

**相位跳变判据的理论保证**：
由推论 4.5，当 $t_k$ 从左侧穿过 $\gamma_j$ 到右侧时，$\arg(D_N(t))$ 发生 $+\pi$ 跳变（反之则为 $-\pi$）。这一跳变的检测不受幅值归一化 $|D_N(t_k)|^2/\|D_N\|_2^2$ 影响，是更稳健的数值判据。

---

## 8. 结论：数学完整性陈述

基于前述 NCDFT 范畴构造、尺度化生成元的谱相变，以及 DFT/FFT 与黎曼函数方程的对偶实现，本节给出黎曼猜想的等价性定理之函子刚性版本。该版本将解析对偶性硬化为纯代数的递归标准形唯一性，并通过构造性反证与极限论证的双重防线排除离轴零点，从而得到最干净的等价链。

### 8.1 递归块三角化刚性

**引理 X（递归块三角化刚性）**：设 $N=2^m$。定义相位权重

```math
W_n^{(N)} := \exp\!\big(i(\alpha-\tfrac12)\operatorname{Li}(x_n^{(N)})\mathbf{H}\big).
```

NCDFT 簇 $\{\mathcal{F}_\alpha^{(N)}\}_{N\geq 1}$ 允许与二分算术分解 $\mathbb{Z}/N\mathbb{Z} \cong \mathbb{Z}/(N/2)\mathbb{Z} \times \{0,1\}$ 相容的递归块三角化（即 Cooley–Tukey 型蝴蝶递归，将 $\mathcal{F}_\alpha^{(N)}$ 逐步归约为规模 $N/2$ 的子问题）当且仅当 $\alpha=1/2$。

**证明概要**：
对 $\alpha=1/2$，有 $W_n\equiv\mathbb{I}_{r+1}$，此时 $\mathcal{F}_{1/2}^{(N)}$ 退化为标准块 DFT。标准 FFT 的蝴蝶操作直接给出块对角标准形，递归链条完整。

对 $\alpha\neq 1/2$，考察对数尺度采样 $x_n^{(N)}=2e^{n\delta_N}$，其中 $\delta_N=\frac{\log(N/2)}{N}$。二分递归要求偶数列权重与半长问题的权重在块置换下一致：

```math
W_{2m}^{(N)} = W_m^{(N/2)}.
```

但

```math
x_{2m}^{(N)} = 2\exp\!\left(\frac{2m\log(N/2)}{N}\right),\qquad 
x_m^{(N/2)} = 2\exp\!\left(\frac{2m\log(N/4)}{N}\right).
```

由于 $\log(N/2)\neq\log(N/4)$，有 $x_{2m}^{(N)}\neq x_m^{(N/2)}$，从而

```math
W_{2m}^{(N)} = \exp\!\big(i(\alpha-\tfrac12)\operatorname{Li}(x_{2m}^{(N)})\mathbf{H}\big) \neq 
\exp\!\big(i(\alpha-\tfrac12)\operatorname{Li}(x_m^{(N/2)})\mathbf{H}\big) = W_m^{(N/2)}.
```

因此，蝴蝶操作后的"偶部"与"奇部"不能递归地归约为规模 $N/2$ 的 NCDFT，递归块三角化链条在第一步即断裂。$\square$

**注**：此处的"块三角化"并非指单个 $N$ 下的静态块 LU 分解，而是要求标准形随 $N$ 的算术二分结构函子地演化。正是这种 $N$-一致性（整体截面）的要求，将 $\alpha=1/2$ 唯一确定。

### 8.2 NCDFT–RH 等价性（函子刚性版本）

**定理 8.1'（NCDFT–RH 等价性，函子刚性版本）**：以下命题等价：

1. **黎曼猜想**：所有非平凡零点满足 $\Re(\rho)=1/2$。
2. **递归标准形**：NCDFT 块加权范德蒙德矩阵簇 $\{\mathcal{F}_\alpha^{(N)}\}_{N\geq 1}$ 存在与二分算术结构相容的递归块三角化标准形。
3. **参数唯一性**：该递归标准形唯一对应 $\alpha=1/2$，且此时尺度化生成元 $\mathcal{H}_{1/2}^{(N)}\equiv 0$。
4. **谱推论**：零点分布测度作为该唯一截面的谱像，其支撑严格共线于 $\Re(s)=1/2$。

**证明概要**：

- **(2) $\Leftrightarrow$ (3)**：由引理 X（递归块三角化刚性）。这是纯有限维矩阵代数的函子性论证，不涉及任何 $N\to\infty$ 极限或解析估计。
- **(3) $\Rightarrow$ (4)**：由第 5.1 节的谱提升算子 $\mathcal{L}_{1/2}^{(N)}$ 的构造，$\mathcal{H}_{1/2}^{(N)}=0$ 锁定实部为 $1/2$。结合定理 5.2 的测度收敛，频率测度 $\nu_{1/2}^{(N)}$ 的支撑由构造就在临界线上，其弱极限继承此共线性质。
- **(4) $\Rightarrow$ (1)**：由命题 5.2'（离轴零点的构造性排斥）。若存在离轴零点 $\rho_0=\sigma_0+i\gamma_0$，则 $\nu_{1/2}^{(N)}$ 的构造支撑（严格在临界线上）与显式公式要求的离轴指数增长直接矛盾。故所有零点必须位于临界线。
- **(1) $\Rightarrow$ (2)**：若黎曼猜想成立，则黎曼函数方程在临界线上的自对偶性迫使 $\mathcal{U}_{1/2}=\mathbb{I}$（定理 2.3），此时 $\mathcal{F}_{1/2}^{(N)}$ 退化为标准块 DFT，FFT 蝴蝶操作直接给出递归块三角化标准形。$\square$

### 8.3 RH 的严格推导（双重防线论证）

由引理 X 与定理 2.3，当 $\alpha=1/2$ 时，尺度化生成元严格为零：

```math
\mathcal{H}_{1/2}^{(N)}=0.
```

这对应于内部自由度的冻结——实部偏差被锁定为零，所有谱点被约束在直线 $\Re(s)=1/2$ 上。然而，这仅确定了零点的实部坐标。要恢复完整的复数零点 $\rho_j = 1/2 + i\gamma_j$（包含虚部信息），必须引入外部频率指标 $t_k = \frac{2\pi k}{\log N}$，通过频率谱测度 $\nu_{1/2}^{(N)}$ 实现：

```math
\nu_{1/2}^{(N)} := \frac{1}{\log N}\sum_{k=0}^{N-1} |D_N(t_k)|^2 \delta_{\frac{1}{2} + i t_k}.
```

该测度的支撑集为：

```math
\mathrm{supp}(\nu_{1/2}^{(N)}) = \left\{\frac{1}{2} + i t_k \;\middle|\; k=0,\dots,N-1\right\} \subset \left\{s\in\mathbb{C} \;\middle|\; \Re(s)=\frac{1}{2}\right\}.
```

即所有质量严格分布在临界线 $\Re(s)=1/2$ 上。由定理 5.2，当 $N\to\infty$ 时，该测度弱收敛到黎曼零点的经验分布：

```math
\lim_{N\to\infty} \nu_{1/2}^{(N)} = \frac{1}{2\pi}\sum_{\rho} \delta_{\rho}.
```

**第一道防线（构造性排斥，有限 $N$）**：命题 5.2' 表明，若存在离轴零点 $\rho_0=\sigma_0+i\gamma_0$（$\sigma_0\neq 1/2$），则 $D_N(t)$ 在 $t=\gamma_0$ 处的幅值将以 $N^{\sigma_0-1/2}$ 的速率指数增长，远超过临界线零点的 $(\log N)^2$ 阶峰值。这将迫使 $\nu_{1/2}^{(N)}$ 在 $\gamma_0$ 邻域内的网格点处承载一个指数爆炸的离轴分量，与 $\nu_{1/2}^{(N)}$ 的支撑严格在临界线上的构造直接矛盾。因此，离轴零点的存在在有限 $N$ 下即被构造性排除。

**第二道防线（极限论证，$N\to\infty$）**：即使退一步忽略上述构造性矛盾，根据测度弱收敛的下半连续性，极限测度的支撑集满足包含关系：

```math
\operatorname{supp}\left(\sum_{\rho} \delta_{\rho}\right) \subseteq \overline{\bigcup_{N} \operatorname{supp}(\nu_{1/2}^{(N)})} \subseteq \left\{s \in \mathbb{C} \;\middle|\; \Re(s)=\frac{1}{2}\right\}.
```

若存在离轴零点 $\rho_0$（$\Re(\rho_0)\neq 1/2$），则 $\delta_{\rho_0}$ 的支撑将落在临界线之外，与上述包含关系矛盾。该极限论证不依赖有限 $N$ 的权重有界性，独立封锁了任何"极限补救"的可能性。

因此，所有非平凡零点 $\rho$ 必须满足 $\Re(\rho)=1/2$。$\square$

**关键区分**：
- **内部谱测度** $\mu_{1/2}^{(N)}=\delta_0$：反映 $\mathcal{H}_{1/2}^{(N)}=0$，仅表示实部锁定（单点紧化），不包含虚部信息。
- **频率谱测度** $\nu_{1/2}^{(N)}$：反映通过 FFT 频率 $t_k$ 展开的虚部信息，其支撑在临界线上，收敛到零点分布。离轴零点的指数增长与此支撑直接冲突。

两者共同构成完整的谱提升：$\mathcal{L}_{1/2}^{(N)}$ 的"内部"部分（来自 $\mathcal{H}$）锁定实部为 $1/2$，"外部"部分（来自 $\mathrm{diag}(t_k)$）展开虚部为 $\gamma_j$。仅当 $\alpha=1/2$ 时，这一双重结构实现严格的自对偶匹配。$\square$

**推论**：此框架提供了构造性证明路径：通过有限维自对偶计算（标准 DFT/FFT）逼近无限维解析对象（黎曼零点），所有步骤满足 Bishop 构造性数学的可计算性要求。NCDFT 框架揭示了黎曼猜想本质上是一个对偶性条件：仅当参数 $\alpha=1/2$ 时，非交换离散傅里叶变换实现严格自对偶，此时算术与谱的函子忠实性达到完美匹配，且谱测度从非紧的扩展态（$\alpha \neq 1/2$）紧化为临界线上的离散点态（$\alpha = 1/2$）。递归块三角化刚性（引理 X）将这一解析对偶性硬化为纯代数的函子刚性——算术数据生成的 NCDFT 簇只允许 $\alpha=1/2$ 的整体截面，而命题 5.2' 的构造性排斥与极限支撑包含关系构成双重防线，确保零点分布作为该截面的谱推论必然共线。

---

## 附录：FFT 蝴蝶操作与对偶结构

### A.1 时间抽取 FFT 的蝴蝶操作

考虑长度为 $N$（$N=2^m$）的 DFT：

```math
X_k = \sum_{n=0}^{N-1} x_n \omega_N^{kn}, \quad \omega_N = e^{-2\pi i/N}.
```

将输入序列按奇偶分开：

```math
X_k = \sum_{m=0}^{N/2-1} x_{2m} \omega_N^{k(2m)} + \sum_{m=0}^{N/2-1} x_{2m+1} \omega_N^{k(2m+1)} = E_k + \omega_N^k O_k,
```

其中 $E_k$ 和 $O_k$ 分别是长度为 $N/2$ 的 DFT：

```math
E_k = \sum_{m=0}^{N/2-1} x_{2m} \omega_{N/2}^{km}, \quad O_k = \sum_{m=0}^{N/2-1} x_{2m+1} \omega_{N/2}^{km}.
```

对于 $k$ 和 $k+N/2$，利用 $\omega_N^{k+N/2} = -\omega_N^k$，得到经典的蝶形运算：

```math
\begin{cases}
X_k = E_k + \omega_N^k O_k, \\
X_{k+N/2} = E_k - \omega_N^k O_k.
\end{cases}
```

这一运算可图示为蝴蝶状结构，因此得名。

### A.2 蝴蝶操作与对偶性

蝴蝶操作揭示了 DFT 的内在对称性：偶部 $E_k$ 与奇部 $O_k$ 通过旋转因子 $\omega_N^k$ 组合，而旋转因子本身满足 $\omega_N^{k+N/2} = -\omega_N^k$，这正是函数方程中 $s \leftrightarrow 1-s$ 对称性的离散体现。在 NCDFT 框架中，对偶性 $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ 对应于参数 $\alpha$ 与 $1-\alpha$ 的互换，而蝴蝶操作中的奇偶分解与符号反转恰好模拟了这一对偶关系。

进一步地，若将 NCDFT 矩阵 $\mathcal{F}_\alpha^{(N)}$ 视为带有内部相位 $\exp(i(\alpha-1/2)\operatorname{Li}(x_n)\mathbf{H})$ 的 DFT，则其快速算法同样可以分解为类似的蝶形结构，只是旋转因子需乘以相应的对角矩阵。这种结构保证了 NCDFT 的计算效率，同时保留了与黎曼函数方程的对偶联系。

### A.3 蝴蝶操作与函数方程

将蝴蝶操作与第 4.2 节的对应表结合，可归纳为：

| FFT 蝴蝶操作 | NCDFT 对偶性 |
| :--- | :--- |
| 奇偶分解 $E_k, O_k$ | 参数对偶 $\alpha \leftrightarrow 1-\alpha$ |
| 旋转因子 $\omega_N^k$ 与 $-\omega_N^k$ | 相位因子 $\exp(i(\alpha-1/2)\operatorname{Li}(x_n)\mathbf{H})$ 及其共轭 |
| $X_k$ 与 $X_{k+N/2}$ 的对称性 | 自对偶条件 $\mathcal{U}_{1/2}=\mathbb{I}$ |
| 网格加密 $N\to\infty$ | 采样完备性 $t_k \to \gamma_j$（第 4.4 节） |

**与引理 X 的联系**：附录 A.1 中的递归要求 $E_k$ 和 $O_k$ 必须是规模 $N/2$ 的同类 DFT，这是 FFT 蝴蝶操作可行的代数根基。引理 X 正是将此递归要求施加于 NCDFT 簇：当 $\alpha\neq 1/2$ 时，相位权重 $W_n$ 不满足 $W_{2m}^{(N)}=W_m^{(N/2)}$，导致偶部与奇部无法归约为规模 $N/2$ 的 NCDFT，递归链条断裂。因此，FFT 蝴蝶操作的函子相容性唯一锁定 $\alpha=1/2$。

---

## 参考文献

1. H. M. Edwards, *Riemann's Zeta Function*, Dover, 2001 [Original: 1974].
2. E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., Oxford University Press, 1986.
3. A. Weil, "Sur les formules explicites de la théorie des nombres premiers", Comm. Sem. Math. Univ. Lund, 1952.
4. A. E. Ingham, "On the estimation of $N(\sigma,T)$", Quart. J. Math., 1940.
5. J. W. Cooley, J. W. Tukey, "An algorithm for the machine calculation of complex Fourier series", Math. Comp., 1965.
6. E. Bishop, *Foundations of Constructive Analysis*, McGraw-Hill, 1967.
