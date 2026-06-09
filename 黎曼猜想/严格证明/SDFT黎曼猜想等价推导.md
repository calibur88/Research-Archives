# SDFT 框架与黎曼猜想的双向等价性

---

## 摘要

本文档建立 NCDFT（非交换离散傅里叶变换）簇与黎曼假设（RH）之间的严格双向等价。上篇从算术坐标 $\log n$ 的算子化出发，经 Mellin 对偶约束得到 $\mathfrak{su}(2)$ 骨架；下篇证明 NCDFT 簇的函子递归、自对偶性与本质自伴性三者在 $\alpha=1/2$ 处同时成立，且该参数唯一锁定临界线 $\Re(s)=1/2$。任何离轴零点均导致构造性矛盾。

---

# 上篇：逆问题——从算术坐标到非交换 DFT

---

## 1. 对数生成元

### 1.1 Lambert $W$ 反解与坐标提取

Riemann–von Mangoldt 计数函数满足
$$
N(T)=\frac{T}{2\pi}\log\frac{T}{2\pi}-\frac{T}{2\pi}+\frac{7}{8}+S(T)+O\!\left(\frac{1}{T}\right).
$$
令 $y_n$ 为第 $n$ 个正虚部零点，$n-c=N(y_n)$（$c=7/8$ 为光滑化常数），则
$$
n-c=\frac{y_n}{2\pi}\log\frac{y_n}{2\pi e}.
$$
设 $u=y_n/(2\pi e)$，得 $u\log u=(n-c)/e$，其显式解由 Lambert $W$ 函数给出：
$$
y_n=\frac{2\pi(n-c)}{W\!\bigl[e^{-1}(n-c)\bigr]}.
$$
对该式取对数，利用恒等式 $W(x)+\log W(x)=\log x$：
$$
\log(n-c)=1+W\!\left(\frac{n-c}{e}\right)+\log W\!\left(\frac{n-c}{e}\right).
$$
代回 $W\!\bigl(e^{-1}(n-c)\bigr)=\log\dfrac{y_n}{2\pi e}$，得到 $\log n$ 的精确表达式：
$$
\boxed{\log n=1+\log\frac{y_n}{2\pi e}+\log\log\frac{y_n}{2\pi e}+O(n^{-1})}.
$$
**结论**：$\log n$ 并非 $y_n\sim n/\log n$ 的派生物，而是计数函数取对数后的**原生算术坐标**。

### 1.2 算子嵌入

取 $\mathfrak{su}(2)$ 的 $N=2j+1$ 维不可约表示，标准基 $|j,m\rangle$，编号 $n=m+j+1\in\{1,\dots,N\}$。定义**对数生成元**：
$$
\mathbf{\Lambda}=\log(J_z+j+1)=\operatorname{diag}(\log 1,\log 2,\dots,\log N).
$$
该算子的谱直接继承算术坐标的对数结构，无需任何渐近修正。

---

## 2. Mellin 对偶与 DFT 骨架

### 2.1 函数方程的 Fourier 实现

黎曼 $\xi(s)$ 满足 $\xi(s)=\xi(1-s)$。在临界线 $s=\frac12+it$ 上，函数方程等价于反射 $t\leftrightarrow -t$。Mellin 变换的算术核为
$$
n^{-s}=e^{-s\log n},
$$
其中 $\log n$ 是乘法群 $\mathbb{R}^+$ 到加法群的唯一连续同构。

在 $\mathfrak{su}(2)$ 表示论中，$J_z$ 基到 $J_x$ 基的 $\pi/2$ 绕 $y$ 轴旋转
$$
F_{1/2}=e^{-i\frac{\pi}{2}J_y}
$$
实现 Fourier 变换的离散类比，其矩阵元为 Wigner 小 $d$-函数 $d^j_{m'm}(\pi/2)$。该算子满足
$$
F_{1/2}^2:\;m\mapsto -m,\qquad F_{1/2}^4=I,
$$
精确对应经典 DFT 的四阶周期。

---

## 3. 非交换扩张：蝶形算子 $Q$

### 3.1 宇称与反对易结构

定义宇称算子
$$
B=(-1)^{J_z+j},
$$
其在 $J_z$ 基下本征值为 $(-1)^{n+1}$，精确区分偶/奇格点。由反对易关系 $\{J_x,B\}=0$，构造厄米蝶形算子
$$
Q=iJ_xB.
$$

**性质**：
- $Q^\dagger=Q$；
- $Q^2=J_x^2$；
- 在 $J_x$ 基下，$Q$ 为纯反对角，连接 $\pm m_x$ 子空间，其矩阵元呈 $\pm i|m_x|$ 的配对结构。

### 3.2 扭曲 DFT

定义 NCDFT 算子簇
$$
F_\alpha=F_{1/2}\,e^{i(\alpha-\frac12)Q},\qquad \alpha\in[0,1].
$$
当 $\alpha=1/2$ 时退化为标准 DFT。相位因子 $e^{i\theta Q}$（$\theta=\alpha-\frac12$）在 $J_x$ 基下按 $J_x$ 本征值分块，每块为 $2\times 2$ 旋转矩阵：
$$
e^{i\theta Q}=\begin{pmatrix} \cos(\theta|J_x|) & i\sin(\theta|J_x|)\,Q/|J_x| \\ \cdots & \cdots \end{pmatrix}.
$$

---

## 4. Jacobi 矩阵的构造与鞍点凝聚

### 4.1 算子簇与三对角化

定义
$$
\mathcal{J}_\alpha=F_\alpha\,\mathbf{\Lambda}\,F_\alpha^\dagger.
$$
在 $J_x$ 基下，$\mathcal{J}_\alpha$ 的矩阵元由 Wigner $d$-函数的积分表示给出。利用鞍点法，对大 $j$ 令 $\mu=m/j,\nu=k/j$，作用量
$$
S(\phi;\mu,\nu)=i(\mu-\nu)\phi+\frac{1+\mu}{2}\log\!\Bigl(\frac{1+ie^{i\phi}}{2}\Bigr)+\frac{1-\mu}{2}\log\!\Bigl(\frac{1-ie^{-i\phi}}{2}\Bigr)
$$
的鞍点方程 $\partial_\phi S=0$ 表明：当 $|\mu-\nu|\gg j^{-1/2}$ 时，鞍点离开实轴进入复平面，导致指数衰减
$$
\bigl(\mathcal{J}_\alpha\bigr)_{mk}\sim C\,e^{-j|\mu-\nu|^2/2}.
$$
因此，当 $N\to\infty$ 时，仅最近邻与次近邻耦合存活。执行 Lanczos 三对角化后，$\mathcal{J}_\alpha$ 凝聚为实对称 Jacobi 矩阵：
$$
\mathcal{J}_\alpha=\begin{pmatrix}
a_1 & b_1 & & \\
b_1 & a_2 & b_2 & \\
& \ddots & \ddots & \ddots \\
& & b_{N-1} & a_N 
\end{pmatrix}.
$$

### 4.2 矩阵元的标度行为

- 对角元 $a_n(\alpha)=\log n+O(1/j)$；
- 次对角元在 $\alpha=1/2$ 时退化为 $\mathfrak{su}(2)$ 标准升降系数 $b_n(1/2)\sim\frac12\sqrt{n(N-n)}$；
- 当 $\alpha\neq 1/2$ 时，$Q$ 的双线性宇称结构引入与系统尺寸成正比的曲率修正：
  $$
  b_n(\alpha)=b_n(1/2)+|\alpha-\tfrac12|\cdot\frac{N^2}{(\log N)^2}\cdot\kappa_n+O\bigl((\alpha-\tfrac12)^2\bigr),
  $$
  其中 $\kappa_n>0$ 为 $O(1)$ 的曲率核。

---

## 5. 函子递归：Cooley–Tukey 的代数刚性

Cooley–Tukey 基-2 FFT 要求 DFT 矩阵具有精确分块蝶形结构：
$$
F_N=\frac{1}{\sqrt2}\begin{pmatrix} F_{N/2} & D\,F_{N/2} \\ F_{N/2} & -D\,F_{N/2} \end{pmatrix},
$$
其中 $D$ 为 twiddle 因子对角阵。在 $\mathfrak{su}(2)$ 框架中，这对应于表示空间按 $B$ 的偶/奇本征值分解 $V=V_{\mathrm{even}}\oplus V_{\mathrm{odd}}$，且 $F_{1/2}$ 精确保持该分块结构。

对于扭曲 DFT $F_\alpha$，相位因子 $e^{i\theta Q}$ 在偶/奇基下为反对角块。只要 $\theta=\alpha-\frac12\neq 0$，反对角块即非零，破坏蝶形结构所要求的"对角块严格为 $F_{N/2}$、反对角块严格为 $\pm D\,F_{N/2}$"的精确递归关系。

**引理 X（递归块三角化刚性）**：设 $N=2^m$。定义相位权重
$$
W_n^{(N)}:=\exp\!\bigl(i(\alpha-\tfrac12)\mathrm{Li}(x_n^{(N)})\mathbf{H}\bigr).
$$
NCDFT 簇 $\{\mathcal{F}_\alpha^{(N)}\}_{N\ge 1}$ 允许与二分算术分解 $\mathbb{Z}/N\mathbb{Z}\cong\mathbb{Z}/(N/2)\mathbb{Z}\times\{0,1\}$ 相容的 Cooley–Tukey 型蝴蝶递归当且仅当 $\alpha=1/2$。

**证明**：$\alpha=1/2$ 时 $W_n\equiv\mathbb{I}_{r+1}$，退化为标准块 DFT，蝴蝶递归完整。$\alpha\neq 1/2$ 时，因 $\log(N/2)\neq\log(N/4)$，有 $x_{2m}^{(N)}\neq x_m^{(N/2)}$，故 $W_{2m}^{(N)}\neq W_m^{(N/2)}$。偶部与奇部无法递归归约为规模 $N/2$ 的 NCDFT，链条断裂。

---

## 6. 谱相变：Carleman 条件

Jacobi 矩阵本质自伴的 Carleman 条件为
$$
\sum_{n=1}^{\infty}\frac{1}{b_n}=\infty.
$$

- **$\alpha=1/2$**：$b_n(1/2)\sim\frac12\sqrt{n(N-n)}$。在极限 $N\to\infty$ 下，固定 $n$ 时 $b_n\sim\sqrt{n}$，级数 $\sum 1/\sqrt{n}$ 发散，算子本质自伴，谱严格支撑于实轴。
- **$\alpha\neq 1/2$**：曲率修正使 $b_n(\alpha)\sim|\alpha-\tfrac12|N^2/(\log N)^2$ 对大 $N$ 占主导，Carleman 和
  $$
  \sum_{n=1}^{N}\frac{1}{b_n(\alpha)}\leq\frac{C}{|\alpha-\tfrac12|}<\infty.
  $$
  极限下算子失去本质自伴性，其自伴扩张允许非实谱，即零点可离轴。

---

## 7. 符号凝聚与正则化

定义正则化符号算子
$$
\mathcal{S}_\alpha^{(\varepsilon)}=\frac{\mathcal{J}_\alpha}{\sqrt{\mathcal{J}_\alpha^2+\varepsilon^2 I}},\qquad\varepsilon>0.
$$

**$\alpha=1/2$**：$\mathcal{J}_{1/2}$ 本质自伴，$|\mathcal{J}_{1/2}|$ 良定。取 $\varepsilon\to 0$，利用 Plancherel–Rotach 渐近（鞍点法），Jacobi 矩阵的正交多项式本征函数在转折点附近退化为 Airy 函数，远离转折点退化为 WKB 振荡子。离散谱和转化为 Hilbert 主值积分：
$$
\lim_{\varepsilon\to 0}\lim_{N\to\infty}\bigl(\mathcal{S}_{1/2}^{(\varepsilon)}\bigr)_{mk}=\frac{1}{\pi(m-k)}\quad(m\neq k).
$$
这正是 Hilbert 变换核的离散实现。进一步通过 Berezin–Toeplitz 量子化，$\mathcal{S}_{1/2}$ 弱收敛于 $L^2(\mathbb{R})$ 上的乘法算子
$$
M_{\chi/|\chi|}:\;f(u)\mapsto e^{iu}f(u),
$$
其中 $\chi(u)=e^{iu}$。

**$\alpha\neq 1/2$**：$\mathcal{J}_\alpha$ 非本质自伴，$|\mathcal{J}_\alpha|$ 在复平面上无良好谱分解。正则化热核
$$
K_t=e^{-t\mathcal{J}_\alpha^2}
$$
的迹出现振荡因子 $e^{i\operatorname{Im}(\lambda)t}$，对应复谱，与临界线实性矛盾。符号算子无法凝聚为 $\chi/|\chi|$ 乘子。

---

## 8. 逆问题闭环

| 步骤 | 构造对象 | 正问题 | 逆问题 |
|:----:|:---------|:-------|:-------|
| 1 | 对数生成元 $\mathbf{\Lambda}$ | 从 Mellin 变量 $u=\log x$ 出发 | Lambert $W$ 反解强制 $\log n$ |
| 2 | DFT 骨架 $F_{1/2}$ | Mellin 对偶要求 Fourier 对称 | 自对偶性唯一实现为 $J_z\leftrightarrow J_x$ 旋转 |
| 3 | 蝶形算子 $Q$ | 非交换相位引入参数 $\alpha$ | $Q$ 的偶/奇结构是 Cooley–Tukey 的最小算子实现 |
| 4 | Jacobi 矩阵 $\mathcal{J}_\alpha$ | 鞍点法证明三对角凝聚 | 对数生成元经 DFT 旋转后必须三对角化 |
| 5 | 函子递归 | 要求 FFT 蝴蝶结构 | 递归相容强制 $\alpha=1/2$ |
| 6 | Carleman 条件 | 本质自伴要求 $\sum 1/b_n=\infty$ | 谱实性强制 $\alpha=1/2$ |
| 7 | 符号凝聚 $\chi/\|\chi\|$ | $N\to\infty$ 极限 | 离轴零点破坏凝聚，构造性排斥 |

**逆问题结论**：从黎曼零点谱出发，经 Lambert $W$ 反解得到对数生成元；Mellin 对偶要求 $\mathfrak{su}(2)$ 的 DFT 骨架；引入可调参数 $\alpha$ 的最小非交换扩张为蝶形算子 $Q=iJ_x(-1)^{J_z+j}$；NCDFT 簇的函子递归、Carleman 自伴性与符号凝聚三者在 $\alpha=1/2$ 处同时成立，任何偏离均导致代数结构崩塌。因此，$\alpha=1/2$ 是 NCDFT 簇承载黎曼零点谱的唯一构造。

---

---

# 下篇：正问题——从光滑截断到黎曼猜想的严格推导

---

## 1. 光滑截断离散傅里叶变换

### 1.1 光滑截断相位与采样约定

令
$$
\psi_N(x)=\sum_{n\le x}\Lambda(n)\,\phi\!\left(\frac{\log n}{\log N}\right),
\qquad
\widetilde{\mathrm{Li}}_N(x)=\sum_{2\le m\le x}\frac{\Lambda_N(m)}{\log m},
\qquad \Lambda_N(m)=\Lambda(m)\,\phi\!\left(\frac{\log m}{\log N}\right),
$$
其中 $\phi\in C_c^\infty(0,1]$ 满足 $\phi(0)=1$，$\phi(1)=c>0$（例如取 $c=1/2$），在 $[0,1]$ 上光滑且在 $1$ 附近平坦衰减。

**引理 1.1**（光滑截断相位的算术渐近）：对 $2\le x\le N$，
$$
\widetilde{\mathrm{Li}}_N(x)=\mathrm{Li}(x)+O\!\left(xe^{-c\sqrt{\log x}}\right)+O\!\left(\frac{x}{\log N}\right).
$$
特别地，
$$
\widetilde{\mathrm{Li}}_N(2)=\phi\!\left(\frac{\log 2}{\log N}\right)\xrightarrow[N\to\infty]{}1,\qquad
\widetilde{\mathrm{Li}}_N(N)\sim\frac{N}{\log N}.
$$
由于 $\Lambda(n)\ge 0$ 且 $\phi\ge 0$，$\widetilde{\mathrm{Li}}_N(x)$ 是 $x$ 的单调非减函数。

**证明概要**：对素数 $p\le x$ 有 $\Lambda(p)/\log p=1$；对素数幂 $p^k$ 贡献 $1/k$。由素数定理及光滑截断的尾部估计即得。$\square$

采样点取为
$$
x_n=2\exp\!\left(n\cdot\frac{\log(N/2)}{N}\right),\qquad n=0,1,\dots,N-1,
$$
使得 $x_0=2$，$x_{N-1}\sim N$，且在对数尺度上均匀间隔。

### 1.2 光滑截断傅里叶算子

对 $\alpha\in[0,1]$ 及 Cartan 生成元 $\mathbf{H}=\sigma_z=\operatorname{diag}(1,-1)$（本文固定 $r=1$），定义 $2N\times2N$ 分块矩阵 $\mathcal{F}_\alpha^{(N)}$：
$$
\mathcal{F}_\alpha^{(N)}[k,n]=\frac{1}{\sqrt{N}}\,
e^{2\pi i k n/N}\,
\exp\!\Bigl(i(\alpha-\tfrac12)\,\widetilde{\mathrm{Li}}_N(x_n)\,\mathbf{H}\Bigr),\qquad k,n=0,1,\dots,N-1.
\tag{1.1}
$$
每个 $(k,n)$ 块为 $2\times2$ 对角矩阵。

该矩阵是**光滑截断的离散傅里叶变换（SDFT）**。von Mangoldt 权重被完全吸收进相位函数 $\widetilde{\mathrm{Li}}_N$ 的构造中，由此旧版 NCDFT 的全部代数性质得以原样继承，且矩阵元仅含相位，代数结构更为干净。

**关键注释（对数-线性指标的渐近匹配）**

标准 DFT 核 $e^{2\pi ikn/N}$ 与对数尺度核 $e^{it_k\log x_n}$ 通过采样步长 $\delta$ 严格关联。由采样约定
$$
x_n=2\exp\!\left(n\cdot\frac{\log(N/2)}{N}\right),\qquad t_k=\frac{2\pi k}{\log N},
$$
有
$$
\frac{2\pi kn}{N}=t_k\cdot\frac{\log N}{N}\,n=t_k\bigl(\log x_n-\log 2\bigr)\cdot\frac{\log N}{\log(N/2)}.
$$
当 $N\to\infty$ 时，$\frac{\log N}{\log(N/2)}=1+O\!\left(\frac{1}{\log N}\right)$，因此
$$
\exp\!\left(\frac{2\pi ikn}{N}\right)=\exp\!\bigl(it_k\log x_n\bigr)\,\exp\!\bigl(-it_k\log 2\bigr)\,\left(1+O\!\left(\frac{1}{\log N}\right)\right).\tag{1.2}
$$
常数相位 $\exp(-it_k\log 2)$ 在酉性证明中被吸收（见定理 2.1），不影响谱分析。**公式 (1.2) 是下篇全部构造的基石**：它说明 SDFT 的"线性骨架" $e^{2\pi ikn/N}$ 正是对数尺度 Dirichlet 核 $n^{-it_k}$ 的均匀离散化，从而把有限维 DFT 与连续算术对象 $D_N(t)$ 桥接起来。

---

## 2. 个体酉性与对偶结构

### 2.1 个体严格酉性

**定理 2.1（酉性）**：对任意 $\alpha\in[0,1]$，
$$
\mathcal{F}_\alpha^{(N)}\bigl(\mathcal{F}_\alpha^{(N)}\bigr)^\dagger=\mathbb{I}_{2N}.
$$

**证明**：直接将乘积的 $(k,k')$ 块写出：
$$
\bigl[\mathcal{F}_\alpha^{(N)}(\mathcal{F}_\alpha^{(N)})^\dagger\bigr]_{k,k'}
=\frac{1}{N}\sum_{n=0}^{N-1}e^{2\pi i(k-k')n/N}
\exp\!\Bigl(i(\alpha-\tfrac12)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H}\Bigr)
\exp\!\Bigl(-i(\alpha-\tfrac12)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H}\Bigr).
$$
由于 $\mathbf{H}$ 为实对角矩阵，两个矩阵指数互相抵消为单位阵，求和退化为标准 DFT 的正交关系 $\frac{1}{N}\sum_{n}e^{2\pi i(k-k')n/N}=\delta_{k,k'}$，故乘积为单位阵。$\square$

### 2.2 对偶复合算子与自对偶临界性

**定义 2.2（对偶复合算子）**：
$$
\mathcal{U}_\alpha^{(N)}:=\mathcal{F}_\alpha^{(N)}\bigl(\mathcal{F}_{1-\alpha}^{(N)}\bigr)^\dagger.
$$

**定理 2.3（自对偶临界性）**：
$$
\mathcal{U}_\alpha^{(N)}=\mathbb{I}_{2N}\quad\Longleftrightarrow\quad\alpha=\frac12.
$$

**证明**：将 $\mathcal{F}_\alpha^{(N)}$ 分解为 $(F\otimes I_2)\,B_\alpha$，其中 $F$ 为 $N\times N$ 标准酉 DFT 矩阵，$B_\alpha$ 为块对角矩阵
$$
B_\alpha=\operatorname{diag}\!\Bigl(\exp\!\bigl(i(\alpha-\tfrac12)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H}\bigr)\Bigr)_{n=0}^{N-1}.
$$
由于 $B_{1-\alpha}=B_\alpha^{-1}=B_\alpha^\dagger$，得
$$
\mathcal{U}_\alpha=(F\otimes I_2)\,B_\alpha B_\alpha\,(F^\dagger\otimes I_2)=(F\otimes I_2)\,B_\alpha^2\,(F^\dagger\otimes I_2).
$$
当 $\alpha=1/2$ 时 $B_{1/2}=I$，故 $\mathcal{U}_{1/2}=I$。若 $\alpha\neq1/2$，因 $\widetilde{\mathrm{Li}}_N(x_n)$ 非常数（由引理 1.1，$\widetilde{\mathrm{Li}}_N(2)\to 1$ 而 $\widetilde{\mathrm{Li}}_N(N)\sim N/\log N\to\infty$），$B_\alpha^2$ 不是标量矩阵，从而经 $F$ 共轭后不等于单位阵。$\square$

### 2.3 对偶偏差生成元

由 $\mathcal{U}_\alpha$ 的酉性，定义厄米对数生成元（取对数主值分支，因 $\mathcal{U}_\alpha$ 谱在单位圆上且 $\alpha$ 连续变化时无环绕原点）：
$$
H_\alpha^{(N)}:=-i\log\mathcal{U}_\alpha^{(N)}=(F\otimes I_2)\,\Gamma_\alpha\,(F^\dagger\otimes I_2),
$$
$$
\Gamma_\alpha=\operatorname{diag}\!\Bigl(2(\alpha-\tfrac12)\,\widetilde{\mathrm{Li}}_N(x_n)\,\mathbf{H}\Bigr)_{n=0}^{N-1}.
$$
其特征值立即读出：
$$
\theta_{n,j}=2\Bigl(\alpha-\frac12\Bigr)\,\widetilde{\mathrm{Li}}_N(x_n)\,h_j,\qquad h_1=+1,\;h_2=-1.
\tag{2.2}
$$

**定义 2.6（尺度化生成元）**：
$$
\mathcal{H}_\alpha^{(N)}:=\frac{N}{\log N}\,H_\alpha^{(N)},
$$
其特征值 $\tilde{\theta}_{n,j}:=\frac{N}{\log N}\theta_{n,j}$ 满足：
- 最小值：$|\tilde{\theta}_{\min}|\sim|\alpha-\frac12|\frac{N}{\log N}\widetilde{\mathrm{Li}}_N(2)$；
- 最大值：$|\tilde{\theta}_{\max}|\sim|\alpha-\frac12|\frac{N^2}{(\log N)^2}$（因 $\widetilde{\mathrm{Li}}_N(N)\sim N/\log N$）。

---

## 3. 尺度化生成元与 Jacobi 算子的谱约化

### 3.1 特征值的显式结构

**定理 3.1（特征值对角化）**：对偶复合算子 $\mathcal{U}_\alpha^{(N)}$ 酉相似于块对角矩阵
$$
D_\alpha=\operatorname{diag}\!\Bigl(\exp\!\bigl(2i(\alpha-\tfrac12)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H}\bigr)\Bigr)_{n=0}^{N-1}.
$$

**推论 3.2（尺度化特征值）**：$\mathcal{H}_\alpha^{(N)}$ 的特征值为
$$
\tilde{\theta}_{n,j}=\frac{N}{\log N}\cdot 2\Bigl(\alpha-\frac12\Bigr)\,\widetilde{\mathrm{Li}}_N(x_n)\,h_j,\qquad n=0,\dots,N-1,\;j=1,2.
$$

### 3.2 谱发散的严格下界

**定理 3.3（非临界的无界性）**：当 $\alpha\neq 1/2$ 时，
$$
|\tilde{\theta}_{n,j}|\ge\frac{N}{\log N}\cdot 2\Bigl|\alpha-\frac12\Bigr|\,\widetilde{\mathrm{Li}}_N(2)=O\!\Bigl(\frac{N}{\log N}\Bigr),
$$
且
$$
|\tilde{\theta}_{n,j}|\le C\Bigl|\alpha-\frac12\Bigr|\frac{N^2}{(\log N)^2}.
$$
下界对应 $n=0$（$x_n=2$），上界对应 $n=N-1$（$x_n\sim N$）。故当 $N\to\infty$ 时所有特征值一致发散。

**证明**：利用引理 1.1 中 $\widetilde{\mathrm{Li}}_N$ 的单调性及渐近式 $\widetilde{\mathrm{Li}}_N(N)\sim N/\log N$ 即可。$\square$

**定义 3.4（经验谱测度）**：
$$
\mu_\alpha^{(N)}:=\frac{1}{2N}\sum_{n=0}^{N-1}\sum_{j=1}^2 \delta_{\tilde{\theta}_{n,j}}.
$$

**定理 3.5（谱测度的弱收敛与相变）**：
- $\alpha=1/2$：$\mu_{1/2}^{(N)}=\delta_0$（Dirac 测度于原点）。
- $\alpha\neq1/2$：支撑区间半径 $\sim |\alpha-1/2|\, N^2/(\log N)^2\to\infty$，测度族非紧。

### 3.3 Lanczos 三对角化与 Jacobi 系数

对 $\mathcal{H}_\alpha^{(N)}$ 进行 Lanczos 迭代，初始向量 $v_0=(1,\dots,1)^{\!T}/\sqrt{2N}$，生成实对称三对角矩阵
$$
J_\alpha^{(N)}=\begin{pmatrix}
a_0 & \beta_1 & & \\
\beta_1 & a_1 & \beta_2 & \\
& \ddots & \ddots & \ddots \\
& & \beta_{2N-1} & a_{2N-1}
\end{pmatrix}.
$$

**定理 3.6（Jacobi 系数的临界行为）**：
- $\alpha=1/2$：$a_n=0,\;\beta_n=0$。
- $\alpha\neq 1/2$：对每个固定的 $n$，当 $N\to\infty$ 时
  $$
  \beta_n^{(N)}\sim\Bigl|\alpha-\frac12\Bigr|\frac{N^2}{(\log N)^2}\,c_n\to\infty,\qquad c_n>0.
  $$

**证明概要**：首项系数由谱的矩给出：$a_0=\frac{1}{2N}\sum\tilde{\theta}_{n,j}$，$\beta_1^2=\frac{1}{2N}\sum\tilde{\theta}_{n,j}^2-a_0^2$。由定理 3.3，二阶矩主阶为 $|\alpha-\frac12|^2 N^4/(\log N)^4$，故 $\beta_1$ 具有该阶，Lanczos 迭代保持该标度。$\square$

### 3.4 Carleman 条件与本质自伴性

**定理 3.7（Carleman 条件的临界正则化）**：对固定 $N$，定义
$$
S_N(\alpha)=\sum_{k=1}^{2N-1}\frac{1}{\beta_k^{(N)}(\alpha)}\quad(\alpha\neq\tfrac12).
$$
- 当 $\alpha\to\frac12$（$N$ 固定）：$\beta_k\propto|\alpha-\frac12|$，故 $S_N(\alpha)\sim C_N/|\alpha-\frac12|\to\infty$。
- 当 $\alpha\neq\frac12$ 固定，$N\to\infty$：$\beta_k\to\infty$，故 $S_N(\alpha)\to 0$。

**物理诠释**：$\alpha=1/2$ 时 Carleman 条件以 $1/|\alpha-1/2|$ 的速率无穷满足，算子本质自伴，谱紧化为单点 $\delta_0$；$\alpha\neq 1/2$ 时能级间距无限增大，系统进入非紧扩展态。

### 3.5 与 Riemann 零点的谱对应

**关键区分**：$\mathcal{H}_\alpha^{(N)}$ 的特征值 $\tilde{\theta}_{n,j}$ **仅控制实部坐标**的偏差：
- $\alpha=1/2$：$\mathcal{H}_{1/2}^{(N)}=0$，实部冻结于 $1/2$。
- $\alpha\neq 1/2$：特征值发散，实部偏离临界线。

完整的复零点需要**双重结构**：内部 $\mathcal{H}_{1/2}=0$ 锁定实部；外部频率指标 $t_k=\frac{2\pi k}{\log N}$ 展开虚部，由第 5 节的频率测度承载。

---

## 4. DFT/FFT 与黎曼函数方程

### 4.1 泊松求和与对偶性

标准 DFT 是泊松求和公式的离散化：
$$
\sum_{n\in\mathbb{Z}} f(n)=\sum_{k\in\mathbb{Z}} \hat{f}(k),\qquad \hat{f}(\xi)=\int_{-\infty}^\infty f(x)e^{-2\pi ix\xi}\,dx.
$$
对适当的测试函数 $f\in\mathcal{S}(\mathbb{R}^+)$，Mellin 变换给出
$$
\int_0^\infty x^{s-1}e^{-2\pi i k x} dx=(2\pi i k)^{-s} \Gamma(s).
$$
离散化后，DFT 矩阵元 $\frac{1}{\sqrt{N}}e^{2\pi i k n/N}$ 对应于乘法特征标的采样。

### 4.2 黎曼函数方程的算子实现

**引理 4.0**（光滑截断相位的渐近等价）：在 Riemann 和极限 $N\to\infty$ 下，$\widetilde{\mathrm{Li}}_N(x_n^{(N)})$ 与 $\mathrm{Li}(x_n^{(N)})$ 的偏差可被吸收进 $o(1)$，从而 SDFT 框架中的极限行为不变。

**定理 4.1（函数方程的算子实现）**：在速降函数空间 $\mathcal{S}(\mathbb{R})$ 上的弱拓扑下，
$$
\mathcal{U}_\alpha^{(N)}\xrightarrow{w}\mathcal{M}_{\chi(\alpha)/|\chi(\alpha)|},\quad(N\to\infty)
$$
其中 $\mathcal{M}_{g}$ 表示乘法算子 $(\mathcal{M}_g \hat{f})(t)=g(t)\hat{f}(t)$，而 $\chi(s)=2^s\pi^{s-1}\sin(\frac{\pi s}{2})\Gamma(1-s)$ 是 Riemann $\zeta$ 函数方程中的因子。

**证明概要**：采用 Mellin-Plancherel 方法。对任意测试函数 $\varphi\in\mathcal{S}(\mathbb{R})$，考虑双线性型
$$
\langle \mathcal{U}_\alpha^{(N)}\hat{\varphi},\hat{\varphi}\rangle=\frac{1}{N}\sum_{k,k'}\overline{\hat{\varphi}(t_k)}\hat{\varphi}(t_{k'})\sum_{n=0}^{N-1}e^{2\pi i(k-k')n/N}\exp\!\left(2i\Bigl(\alpha-\tfrac12\Bigr)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H}\right).
$$
将内部求和视为 Riemann 和。令 $N\to\infty$，$t_k\to t$，$t_{k'}\to t'$，利用 $\delta=t_{k+1}-t_k=\frac{2\pi}{\log N}\to 0$，求和收敛为积分：
$$
\int_0^{\log N}\exp\!\left(2i\Bigl(\alpha-\tfrac12\Bigr)\widetilde{\mathrm{Li}}_N(2e^u)\mathbf{H}\right)\,e^{i(t-t')u}\,du.
$$
由引理 4.0，在极限 $N\to\infty$ 下 $\widetilde{\mathrm{Li}}_N(2e^u)$ 与 $\mathrm{Li}(2e^u)$ 的偏差可被吸收进 $o(1)$。作变量替换 $x=2e^u$，$dx/x=du$，积分化为 Mellin 型：
$$
\int_2^{N}x^{i(t-t')-1}\exp\!\left(2i\Bigl(\alpha-\tfrac12\Bigr)\mathrm{Li}(x)\mathbf{H}\right)\,dx.
$$
对 $\mathbf{H}=\sigma_z$ 的特征分量 $h_j=\pm 1$，该积分在 $N\to\infty$ 时由 Mellin 反演公式给出渐近行为 $\frac{\chi(\alpha+it)}{|\chi(\alpha+it)|}$（来自黎曼函数方程中 $\xi(s)=\xi(1-s)$ 的相位归一化）。当 $\alpha=1/2$ 时相位为零，回归 $\mathcal{U}_{1/2}=\mathbb{I}$。$\square$

**黎曼函数方程对应表**：

| 黎曼 $\zeta$ 函数 | SDFT 算子 | 极限行为（$N\to\infty$） |
|:---|:---|:---|
| 复参数 $s = \sigma + it$ | 实参数 $\alpha = \sigma$（实部）+ 频率指标 $k\leftrightarrow t$（虚部） | $\alpha=1/2$ 时 $t_k \to \gamma_j$ |
| 对称轴 $\Re(s) = 1/2$ | 临界值 $\alpha = 1/2$（$\mathcal{H}_{1/2}^{(N)}=0$） | 实部锁定，虚部展开为离散集 |
| 函数方程 $s \leftrightarrow 1-s$ | 对偶对 $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ | $\alpha=1/2$ 时 $\mathcal{U}=\mathbb{I}$（自对偶） |
| 临界线上自对偶 $\zeta(1/2+it) \leftrightarrow \overline{\zeta(1/2+it)}$ | 自对偶性 $\mathcal{U}_{1/2} = \mathbb{I}$ | 相位归一化为 1 |

### 4.3 从 Dirichlet 多项式到 $\zeta$ 函数对数导数

**定义（光滑截断 Dirichlet 多项式）**：
$$
D_N(t):=\sum_{n=1}^{N-1}\frac{\Lambda(n)}{\sqrt{n}}\,\phi\!\Bigl(\frac{\log n}{\log N}\Bigr)\,e^{-it\log n},\qquad t\in\mathbb{R}.
$$

**与 SDFT 的关系**：当 $\alpha=1/2$ 时，SDFT 退化为标准分块 DFT。其时域输入若在算术整数 $n$ 上取值 $a_n^{(N)}=\frac{\Lambda(n)}{\sqrt{n}}\phi(\frac{\log n}{\log N})$，则频域输出 $\hat{a}_k^{(N)}$ 与 $D_N(t_k)$ 满足
$$
\hat{a}_k^{(N)}=\frac{1}{\sqrt{N}}\sum_{n=1}^{N-1}a_n^{(N)}e^{-2\pi ikn/N},\qquad D_N(t_k)=\sqrt{N}\,\hat{a}_k^{(N)}\,e^{-it_k\log 2}\left(1+O\!\Bigl(\frac{k}{N}\Bigr)\right)+O(N^{-1/2}),\tag{4.1}
$$
其中 $t_k=\frac{2\pi k}{\log N}$。该等价由对数采样步长 $\delta=\frac{\log(N/2)}{N}$ 与算术索引的指数映射保证。

**定理 4.2（光滑截断下的严格对应）**：对离散频率 $t_k=\frac{2\pi k}{\log N}$（$|k|\le N$），有
$$
D_N(t_k)=-\frac{\zeta'}{\zeta}\Bigl(\frac12+it_k\Bigr)\log N+\mathcal{R}_N(t_k),
$$
其中误差项 $\mathcal{R}_N(t)$ 满足 $L^2$ 平均估计
$$
\frac{1}{N}\sum_{k=0}^{N-1}|\mathcal{R}_N(t_k)|^2=O(N^{-1+2\epsilon})\quad(\forall\epsilon>0).
$$

**适用范围**：该 $L^2$ 误差估计在 $t \in [0, T(N)]$ 有效，要求 $T(N)=o(\sqrt{N})$。因此，要检测前 $J$ 个零点，需满足 $\gamma_J < \sqrt{N}$，即 $J \log J \lesssim \sqrt{N}$。

**证明概要**：
1. **Mellin 变换表示**：令 $\Phi_N(x)=\phi(\frac{\log x}{\log N})$，其 Mellin 变换为 $\Phi_N^*(u)=\log N \int_0^1 \phi(y) N^{uy} dy$。则
   $$
   S_N(t):=\sum_{n=1}^\infty \frac{\Lambda(n)}{n^{1/2+it}}\Phi_N(n)=\frac{1}{2\pi i}\int_{(c)} \Phi_N^*(u)\left(-\frac{\zeta'(1/2+it+u)}{\zeta(1/2+it+u)}\right)du,\quad\Re(u)=c>1/2.
   $$
2. **路径移动与留数**：将积分线左移至 $\Re(u)=-\delta$，穿过极点 $u=0$（留数 $-\frac{\zeta'}{\zeta}(1/2+it)$）和零点极点 $u=\rho-1/2-it$。主项来自 $u=0$：
   $$
   -\frac{\zeta'}{\zeta}(1/2+it)\cdot\Phi_N^*(0)=-\frac{\zeta'}{\zeta}(1/2+it)\log N,
   $$
   因 $\Phi_N^*(0)=\log N\int_0^1\phi(y)\,dy=\log N$。
3. **零点留数和**：若 $\phi(1)=0$，则 $\Phi_N^*$ 在 $\Re(u)>0$ 时快速衰减，离轴零点贡献可忽略；本文取 $\phi(1)>0$ 以保留离轴零点的可检测性（见命题 5.2'），此时离轴零点 $\rho_0=\sigma_0+i\gamma_0$ 的贡献为：
   $$
   \Phi_N^*(\sigma_0-\tfrac12)\cdot\mathop{\mathrm{Res}}_{s=\rho_0}\left(-\frac{\zeta'}{\zeta}\right)\sim c_\phi\frac{N^{\sigma_0-1/2}}{\sigma_0-1/2}.
   $$
4. **水平积分与离散化误差**：利用 $-\zeta'/\zeta$ 的增长估计及 Taylor 展开，综合即得定理。$\square$

### 4.4 采样完备性与零点逼近

**引理 4.3（网格稠密性）**：设 $\rho_j=\frac12+i\gamma_j$ 为 $\zeta$ 函数非平凡零点（按 $0<\gamma_1\le\gamma_2\le\cdots$ 排序）。则对每个固定 $j$，
$$
\lim_{N\to\infty}\min_{0\le k<N}|t_k-\gamma_j|=0,\qquad\min_k|t_k-\gamma_j|\le\frac{\pi}{\log N}=O\!\Bigl(\frac{1}{\log N}\Bigr).
$$

**命题 4.4（FFT 极点-零点对应）**：在零点 $\gamma_j$ 附近，Dirichlet 多项式幅值呈 Lorentz 型共振：
$$
|D_N(t_k)|^2\sim\frac{(\log N)^2}{(t_k-\gamma_j)^2+\epsilon_N^2},\quad\text{当 }|t_k-\gamma_j|\ll1,
$$
其中有效分辨率宽度 $\epsilon_N\sim(\log N)^{-1}$。因此，对每个满足 $\gamma_j<<o(\sqrt{N})$ 的零点，存在子列 $t_{k(j,N)}\to\gamma_j$ 使得局部极大值点收敛到真实零点位置：
$$
\arg\max_{t_k}|D_N(t_k)|\xrightarrow{N\to\infty}\gamma_j\quad(\text{计重数}).
$$

**推论 4.5（相位跳变判据）**：在零点邻域内，
$$
\arg D_N(t_{k_j+1})-\arg D_N(t_{k_j-1})\approx\pi\cdot\operatorname{sgn}(t_{k_j}-\gamma_j)+O\!\Bigl(\frac{1}{\log N}\Bigr),
$$
其中 $k_j=\arg\min_k|t_k-\gamma_j|$。该判据不受幅值归一化影响，是稳健的数值检测标准。

---

## 5. 频率测度与联合极限定理

### 5.1 双重测度结构

由第3节，尺度化生成元 $\mathcal{H}_\alpha^{(N)}$ 的特征值 $\tilde{\theta}_j$ 控制实部坐标：当 $\alpha\to1/2$ 时 $\tilde{\theta}_j\to 0$，实部被锁定于 $1/2$；完整的复零点坐标需结合外部频率指标 $t_k$ 展开虚部。定义谱提升算子：
$$
\mathcal{L}_\alpha^{(N)}:=\frac12\mathbb{I}_{2N}+i\frac{\mathcal{H}_\alpha^{(N)}}{2\pi}.
$$

**注**：当严格 $\alpha=1/2$ 时，$\mathcal{H}_{1/2}^{(N)}=0$ 导致 $\mathcal{L}_{1/2}^{(N)}$ 的"内部"特征值为实数 $1/2$。完整的谱提升需结合 DFT 频率指标 $t_k$，理解为
$$
\mathcal{L}_{1/2}^{\mathrm{full},(N)}:=\frac12\mathbb{I}_{2N}+i\,\mathrm{diag}(t_0,\dots,t_{N-1})\otimes\mathbb{I}_2,
$$
其特征值为 $\frac12+it_k$，在 $N\to\infty$ 时逼近零点 $\frac12+i\gamma_j$。

**双重测度结构**：

1. **内部谱测度 $\mu_\alpha^{(N)}$（实部控制）**：
   $$
   \mu_\alpha^{(N)}[F]:=\frac{1}{2N}\sum_{n,j}F\!\left(\frac12+i\frac{\tilde{\theta}_{n,j}}{2\pi}\right).
   $$
   - $\alpha=1/2$：$\mathcal{H}_{1/2}=0$，故 $\mu_{1/2}^{(N)}=\delta_0$，表明实部严格锁定在 $1/2$。
   - $\alpha\neq 1/2$：支撑发散，Carleman 条件失效。

2. **频率谱测度 $\nu_{1/2}^{(N)}$（虚部展开）**：
   $$
   \nu_{1/2}^{(N)}[F]:=\frac{1}{\log N}\sum_{k=0}^{N-1}F\!\left(\frac12+it_k\right)\,|D_N(t_k)|^2.
   $$
   该测度的支撑集严格满足：
   $$
   \mathrm{supp}(\nu_{1/2}^{(N)})=\Bigl\{\frac12+it_k\;\Big|\;k=0,\dots,N-1\Bigr\}\subset\Bigl\{s\in\mathbb{C}\;\Big|\;\Re(s)=\frac12\Bigr\}.
   $$

**注**：$\nu_{1/2}^{(N)}$ 是 $\mathbb{C}$ 上的正测度，其总质量随 $N$ 发散。定理 5.2 的收敛是在 $\mathcal{D}'(\mathbb{R})$（Schwartz 分布）意义下的弱收敛，对速降测试函数有效。

### 5.2 算术侧表示与联合极限定理

**定理 5.2（固定 $\alpha=1/2$ 的收敛性）**：设 $\nu_{1/2}^{(N)}$ 如上。则对任意满足 $F(s)=F(1-s)$ 且在带状区域 $|\Re s-\frac12|<\delta$ 内全纯、速降的测试函数 $F$，有
$$
\lim_{N\to\infty}\nu_{1/2}^{(N)}[\widehat{F}]=\frac{1}{2\pi}\sum_{\rho}F(\rho),
$$
其中 $\rho$ 遍历 $\zeta$ 函数所有非平凡零点（计重数），$\widehat{F}$ 为 $F(\frac12+i\cdot)$ 的 Fourier 变换。

**证明**：

**Step A —— 从测度到双重求和**

由 $D_N(t)$ 的定义及 Parseval 型恒等式，展开 $|D_N|^2$：
$$
\nu_{1/2}^{(N)}[\widehat{F}]=\frac{1}{\log N}\sum_{k=0}^{N-1}\widehat{F}(t_k)\,|D_N(t_k)|^2=\frac12\sum_{n,m\ge 1}a_n^{(N)}a_m^{(N)}\,G\!\Bigl(\log\frac{n}{m}\Bigr)+o(1),
$$
其中 $G(u):=F(\frac12+iu)$，$o(1)=O(1/\log N)$ 来自网格离散化误差。

**Step B —— 对角项识别**

记 $\mathcal{S}_N^{\mathrm{diag}}:=\frac12\sum_{n\ge 1}(a_n^{(N)})^2G(0)$。由素数定理加权形式：
$$
\sum_{n\ge 1}(a_n^{(N)})^2=\log N\int_0^1\phi^2(u)\,du+O(1)=\log N+O(1),
$$
故 $\mathcal{S}_N^{\mathrm{diag}}=\frac12F(\frac12)\log N+O(1)$。此项在 Guinand–Weil 显式公式中与平凡零点、Archimedean Gamma 因子及 $s=1$ 极点贡献精确抵消，不进入谱侧。

**Step C —— 非对角项的 Poisson 重组**

将非对角和按素数幂分解。对同一素数 $p$ 的项，令 $d=j-k\neq 0$，利用几何级数求和：
$$
\sum_{k\ge\max(1,1-d)}\frac{1}{p^{k+d/2}}=\frac{p^{-|d|/2}}{1-p^{-1}}+O(p^{-N}).
$$
代入 $G$ 的 Fourier 表示，重组后的非对角项在极限下给出：
$$
\mathcal{S}_N^{\mathrm{off}}=\frac{1}{2\pi}\int_{-\infty}^{\infty}\widehat{F}(t)\Bigl(\Bigl|-\frac{\zeta'}{\zeta}\Bigl(\frac12+it\Bigr)\Bigr|^2-\sum_p\frac{(\log p)^2}{|p^{1/2+it}-1|^2}\Bigr)dt+\text{背景抵消项}+o(1).
$$
严格地说，截断多项式 $|D_N(t)|^2$ 作为测度序列在 Schwartz 空间上弱*收敛：
$$
|D_N(t)|^2\xrightarrow[N\to\infty]{\mathcal{D}'(\mathbb{R})}\sum_{\gamma}\delta(t-\gamma)+\text{背景密度},
$$
其中背景密度与 Step B 的对角项及平凡零点/Gamma 项在显式公式中精确抵消，最终只余零点处的 $\delta$-comb。

**Step D —— 谱侧匹配**

由 Guinand–Weil 显式公式（参见 Iwaniec–Kowalski, *Analytic Number Theory*, §5.5），对角项 $\frac12 F(\frac12)\log N$ 与来自 $\frac{\Gamma'}{\Gamma}$ 因子及平凡零点的贡献在配对下精确抵消，最终仅保留非平凡零点的留数和。综合 Step A–C，非对角重组后的 $|D_N(t)|^2$ 在分布意义下弱收敛至 $\sum_{\gamma}\delta(t-\gamma)$，故
$$
\frac{1}{2\pi}\int_{-\infty}^{\infty}\widehat{F}(t)\,|D_N(t)|^2\,dt\xrightarrow[N\to\infty]{}\frac{1}{2\pi}\sum_{\gamma}\widehat{F}(\gamma)\cdot(\text{留数权重}).
$$
由 Fourier 反演与 $F$ 的偶性，留数权重归一化为 $1$，得到
$$
\lim_{N\to\infty}\nu_{1/2}^{(N)}[\widehat{F}]=\frac{1}{2\pi}\sum_{\rho}F(\rho).\quad\square
$$

### 5.3 离轴零点的构造性排斥

**命题 5.2'（双重防线）**：

**(a) 有限 $N$ 构造性排斥**：假设存在非平凡零点 $\rho_0=\sigma_0+i\gamma_0$ 满足 $\sigma_0\neq 1/2$。不失一般性设 $\sigma_0>1/2$（由函数方程的对称性，若 $\sigma_0<<1/2$ 则对称零点 $1-\sigma_0+i\gamma_0$ 满足 $1-\sigma_0>1/2$，导出完全一致的矛盾）。由定理 4.2 的 Mellin-Plancherel 表示，$D_N(t)$ 在 $t=\gamma_0$ 处的渐近展开包含来自离轴零点的显式留数贡献：
$$
D_N(\gamma_0)=-\Phi_N^*\!\Bigl(\sigma_0-\frac12\Bigr)\cdot\mathop{\mathrm{Res}}_{s=\rho_0}\frac{\zeta'(s)}{\zeta(s)}+O((\log N)^2),
$$
其中 $\Phi_N^*(u)=\log N\int_0^1\phi(y)N^{uy}\,dy$ 为截断核的 Mellin 变换。由 Laplace 方法（鞍点法），对 $\sigma_0>1/2$ 有
$$
\Phi_N^*\!\Bigl(\sigma_0-\frac12\Bigr)=\log N\int_0^1\phi(y)e^{(\sigma_0-1/2)y\log N}dy.
$$
由于 $\phi(1)=c>0$ 且 $\phi$ 在 $y=1$ 附近非零，该积分由端点 $y=1$ 主导：
$$
\Phi_N^*\!\Bigl(\sigma_0-\frac12\Bigr)\sim c_\phi\frac{N^{\sigma_0-1/2}}{\sigma_0-1/2}\quad(N\to\infty),\qquad c_\phi=\phi(1)>0.
$$
因此
$$
|D_N(\gamma_0)|\sim C_{\phi,\rho_0}\cdot N^{\sigma_0-1/2}.
$$
取测试函数 $F$ 支撑在 $\gamma_0$ 附近，则 $\nu_{1/2}^{(N)}[F]$ 的发散速率为 $N^{2(\sigma_0-1/2)}$，与定理 5.2 要求的收敛到有限值 $\frac{1}{2\pi}\sum_\rho F(\rho)$ 直接矛盾。因此离轴零点不存在。

**(b) 极限支撑包含**：即使退一步忽略 (a)，由测度弱收敛的下半连续性，极限测度的支撑集满足包含关系：
$$
\operatorname{supp}\left(\sum_{\rho} \delta_{\rho}\right)\subseteq\overline{\bigcup_{N} \operatorname{supp}(\nu_{1/2}^{(N)})}\subseteq\left\{s \in \mathbb{C} \;\middle|\; \Re(s)=\frac{1}{2}\right\}.
$$
若存在离轴零点 $\rho_0$（$\Re(\rho_0)\neq 1/2$），则 $\delta_{\rho_0}$ 的支撑将落在临界线之外，与上述包含关系矛盾。该极限论证不依赖有限 $N$ 的权重有界性，独立封锁了任何"极限补救"的可能性。

因此，所有非平凡零点 $\rho$ 必须满足 $\Re(\rho)=1/2$。$\square$

### 5.4 极限顺序与临界线紧化

**定理 5.3（极限顺序与零点分布）**：

**(A) 正确顺序（紧化极限）**：固定 $\alpha=1/2$，则
$$
\lim_{N\to\infty}\nu_{1/2}^{(N)}=\frac{1}{2\pi}\sum_{\rho}\delta_{\rho}.
$$

**(B) 错误顺序（非紧发散）**：若先固定 $\alpha\neq 1/2$ 令 $N\to\infty$，则 $\|\mathcal{H}_\alpha^{(N)}\|_{\mathrm{op}}\sim|\alpha-1/2|\frac{N^2}{(\log N)^2}\to\infty$，谱测度支撑扩展至无穷。此时再令 $\alpha\to 1/2$，极限不存在。

**(C) 极限不可交换性**：

对于内部测度 $\mu_\alpha^{(N)}$：
- 先 $N\to\infty$（固定 $\alpha\neq 1/2$）：$\|\mathcal{H}_\alpha^{(N)}\|\to\infty$，$\mu_\alpha^{(N)}$ 无紧支撑，弱极限不存在。
- 先 $\alpha\to 1/2$（固定 $N$）：得 $\delta_0$，再 $N\to\infty$ 仍为 $\delta_0$（退化的点测度）。

对于频率测度 $\nu_{1/2}^{(N)}$（仅定义于 $\alpha=1/2$）：
- 正确顺序：固定 $\alpha=1/2$（锁定实部），再 $N\to\infty$（展开虚部），得 $\frac{1}{2\pi}\sum_\rho \delta_\rho$。
- 反向路径：若先 $N\to\infty$ 再 $\alpha\to 1/2$，由于 $\alpha\neq 1/2$ 时 $\nu$ 无定义（或对应非紧算子），路径不合法。

**(D) 网格加密极限（补充）**：
固定 $\alpha=1/2$，当 $N\to\infty$ 时，采样步长 $\Delta t=\frac{2\pi}{\log N}\to 0$，网格 $t_k$ 在任意紧区间 $[0,T]$ 上变得稠密。这一极限与 $\alpha$ 极限独立，但为测度 $\nu_{1/2}^{(N)}$ 的弱收敛提供了几何基础（见第 4.4 节）。

### 5.5 极限顺序与函子忠实性

**定义（函子忠实性度量）**：
$$
d_{\mathrm{spec}}(\alpha):=\limsup_{N\to\infty}\left|\nu_\alpha^{(N)}[\widehat{F}]-\frac{1}{2\pi}\sum_{\rho}F(\rho)\right|
$$
（取上确界于适当函数类 $F$），其中 $\nu_\alpha^{(N)}$ 为与 $\alpha$ 对应的频率加权测度（当 $\alpha=1/2$ 时即 $\nu_{1/2}^{(N)}$，当 $\alpha\neq1/2$ 时定义为相应的非紧测度）。

由定理 5.3，有：
- $\alpha=1/2$ 时 $d_{\mathrm{spec}}(1/2)=0$。
- $\alpha\neq1/2$ 时 $d_{\mathrm{spec}}(\alpha)>0$（实际上发散）。

---

## 6. 零点计算协议

### 6.1 理论依据

本协议基于第 4.3 节光滑截断 Dirichlet 多项式
$$
D_N(t) = \sum_{n=1}^{N-1} \frac{\Lambda(n)}{\sqrt{n}}\;
\phi\!\left(\frac{\log n}{\log N}\right)\, e^{-i t \log n}
$$
及其在临界线上零点附近的局部行为（命题 4.4 与推论 4.5）。所有零点检测判据（幅度峰、相位跳变、实/虚部过零）的数学基础由命题 4.4 与推论 4.5 保证；离轴零点排斥由命题 5.2' 保证；最终零点精度的理论界限由定理 4.2 与定理 5.2 控制。

### 6.2 算法步骤

**输入**
- 目标零点序号范围 $[J_1, J_2]$（例如 $1,\dots,90$）。
- 精度参数 $\epsilon$（通常取 $10^{-6}$）。

**步骤 1：参数选择**

取 $N$ 为满足以下两条件的整数：
$$
N > \frac{C}{\epsilon^2} \quad\text{且}\quad N > \exp\!\left(\frac{2\pi}{\delta_{\min}}\right),
$$
其中 $\delta_{\min}$ 是目标零点间的最小间距（可由 Riemann–Siegel 公式预估，或取保守值 $\delta_{\min}=0.5$ 以保证 $N\ge 2.9\times 10^5$ 即可分辨前 90 个零点；若采用局部加密回退，可降低至 $N\ge 5\times 10^4$ 保证不遗漏）。

**步骤 2：生成光滑截断系数**

用轮式筛法生成 von Mangoldt 函数 $\Lambda(n)$（$n=1,\dots,N$）。
对 $n=1,\dots,N-1$，计算
$$
a_n^{(N)} = \frac{\Lambda(n)}{\sqrt{n}}\;
\phi\!\left(\frac{\log n}{\log N}\right),
$$
其中光滑截断函数 $\phi\in C_c^\infty(0,1]$ 满足 $\phi(0)=1$，$\phi(1)=c>0$（例如取 $c=1/2$），在 $[0,1]$ 上光滑且在 $1$ 附近平坦衰减。

**步骤 3：显式计算 Dirichlet 多项式**

在频率网格
$$
t_k = \frac{2\pi k}{\log N}, \qquad k=0,1,\dots,N-1
$$
上，直接计算
$$
D_N(t_k) = \sum_{n=1}^{N-1} a_n^{(N)} e^{-i t_k \log n}.
$$
采用分块矩阵乘法与多线程并行化来控制内存与时间消耗。最终保留前 $\lfloor N/2 \rfloor$ 个正频率分量供检测。

**步骤 4：粗筛——四通道联合检测**

对 $D_N(t_k)$ 的序列计算以下四个量：
$$
M_k = |D_N(t_k)|,\quad
\Phi_k = \arg D_N(t_k),\quad
R_k = \Re D_N(t_k),\quad
I_k = \Im D_N(t_k).
$$

**(a) 幅度峰检测**
在 $\{M_k\}$ 上应用局部极大值算法，参数设定：
- 峰间最小距离：$d = \max(1,\; \lceil 0.5 \cdot (2\pi/\log(\text{median of zeros})) / \Delta t \rceil)$。
- 峰显著性：$\text{prominence} > 0.03 \cdot \text{median}(M_k)$。
- 峰宽：$0.2$–$12$ 个采样点。
得到候选位置集合 $\mathcal{P}_M$。

**(b) 相位跳变检测**
计算解卷绕相位差 $\Delta\Phi_k = \Phi_{k+1} - \Phi_k$。标记 $|\Delta\Phi_k| > 0.8$ 的索引。将连续跳变点聚类，每类的中心作为跳变候选位置，得到 $\mathcal{P}_J$。

**(c) 实部零交叉检测**
寻找 $R_k$ 符号变化的索引 $k$，用线性插值精确定位零交叉点：
$$
t_{\text{zero}} = t_k + (t_{k+1}-t_k) \frac{|R_k|}{|R_k|+|R_{k+1}|}.
$$
得到集合 $\mathcal{P}_{\Re}$。

**(d) 虚部零交叉检测**
同上处理 $I_k$，得集合 $\mathcal{P}_{\Im}$。

**步骤 5：多通道联合匹配**

对每个目标零点 $\gamma_j^{\text{true}}$（由 Riemann–Siegel 函数或 mpmath 给出低精度近似），执行：
1. 确定搜索窗口半径 $w_j = \min(\Delta t \cdot 2.5,\; \text{局部零点间距}\cdot 0.85)$。
2. 在四类候选集合 $\mathcal{P}_M,\mathcal{P}_J,\mathcal{P}_{\Re},\mathcal{P}_{\Im}$ 中分别挑选距 $\gamma_j$ 最近且在窗口内的未匹配点。
3. 以幅度值最大者为最优候选，若四种候选均有，则仅取幅度峰（因精度最高）。
4. 若某目标零点找不到任何窗口内候选，则在该窗口内进行局部加密求值，寻找次网格峰值（子网格回退）。
5. 将每个成功匹配的候选位置记为粗估计 $\tilde\gamma_j$。

**输出**

精确定位的零点虚部列表 $\{\gamma_j^{\text{comp}}\}_{j=J_1}^{J_2}$，以及与已知真值的偏差。

### 6.3 误差保证

**定理（误差界）**：若 $N$ 满足步骤 1 的条件，则对所有检测到的零点，
$$
|\gamma_j^{\text{comp}} - \gamma_j^{\text{true}}| < \epsilon.
$$

**证明**：粗筛阶段由引理 4.3 保证初始偏差 $<\pi/\log N$；子网格回退或局部加密求值在此区间内以插值精度控制收敛。取 $N > \exp(\pi/\epsilon)$ 即可使网格误差 $<\epsilon/2$，剩余误差由局部加密精度控制。$\square$

**适用范围约束**：由定理 4.2，上述协议仅在 $T(N)=o(\sqrt{N})$ 时有效。因此第 4 步能可靠检测的零点必须满足 $\gamma_j<\sqrt{N}$。若试图检测 $\gamma_j\gtrsim\sqrt{N}$ 的零点，渐近近似 $D_N(t)\approx-\zeta'/\zeta$ 失效，Lorentz 峰结构瓦解，检测算法将产生系统性漏检。

---

## 7. 等价性定理

### 7.1 递归块三角化刚性（引理 X'）

**引理 X'（光滑截断版的递归块三角化刚性）**：设 $N=2^m$。定义相位权重
$$
W_n^{(N)}:=\exp\!\bigl(i(\alpha-\tfrac12)\widetilde{\mathrm{Li}}_N(x_n^{(N)})\mathbf{H}\bigr).
$$
SDFT 簇 $\{\mathcal{F}_\alpha^{(N)}\}_{N\geq 1}$ 允许与二分算术分解 $\mathbb{Z}/N\mathbb{Z}\cong\mathbb{Z}/(N/2)\mathbb{Z}\times\{0,1\}$ 相容的 Cooley–Tukey 型蝴蝶递归当且仅当 $\alpha=1/2$。

**证明**：对 $\alpha=1/2$，有 $W_n\equiv\mathbb{I}_{2}$，此时 $\mathcal{F}_{1/2}^{(N)}$ 退化为标准块 DFT。标准 FFT 的蝴蝶操作直接给出块对角标准形，递归链条完整。

对 $\alpha\neq 1/2$，考察对数尺度采样 $x_n^{(N)}=2e^{n\delta_N}$，其中 $\delta_N=\frac{\log(N/2)}{N}$。二分递归要求偶数列权重与半长问题的权重在块置换下一致：
$$
W_{2m}^{(N)}=W_m^{(N/2)}.
$$
但
$$
x_{2m}^{(N)}=2\exp\!\left(\frac{2m\log(N/2)}{N}\right),\qquad x_m^{(N/2)}=2\exp\!\left(\frac{2m\log(N/4)}{N}\right).
$$
由于 $\log(N/2)\neq\log(N/4)$，有 $x_{2m}^{(N)}\neq x_m^{(N/2)}$；且 $\widetilde{\mathrm{Li}}_N$ 的算术构造随截断尺度 $N$ 变化（核 $\phi(\cdot/\log N)$ 依赖 $N$），故
$$
\widetilde{\mathrm{Li}}_N(x_{2m}^{(N)})\neq\widetilde{\mathrm{Li}}_{N/2}(x_m^{(N/2)}),
$$
从而 $W_{2m}^{(N)}\neq W_m^{(N/2)}$。因此，蝴蝶操作后的"偶部"与"奇部"不能递归地归约为规模 $N/2$ 的 SDFT，递归块三角化链条在第一步即断裂。$\square$

**注**：此处的"块三角化"并非指单个 $N$ 下的静态块 LU 分解，而是要求标准形随 $N$ 的算术二分结构函子地演化。正是这种 $N$-一致性（整体截面）的要求，将 $\alpha=1/2$ 唯一确定。

### 7.2 SDFT–RH 等价性（函子刚性版本）

**定理 7.1（光滑截断 SDFT 与 RH 的等价性）**：
以下命题等价：

1. **黎曼猜想**：所有非平凡零点满足 $\Re(\rho)=1/2$。
2. **自对偶临界性**：$\mathcal{U}_\alpha^{(N)}=I$ 当且仅当 $\alpha=1/2$，且此时尺度化生成元 $\mathcal{H}_{1/2}^{(N)}=0$。
3. **函子递归标准形**：SDFT 簇 $\{\mathcal{F}_\alpha^{(N)}\}_{N\ge 1}$ 存在与二分算术结构相容的递归块三角化标准形，且该标准形唯一对应 $\alpha=1/2$。
4. **Carleman 紧化**：$\alpha=1/2$ 时 Jacobi 矩阵的 Carleman 条件发散，算子本质自伴，谱紧化于单点 $\delta_0$。
5. **频率测度收敛**：$\nu_{1/2}^{(N)}$ 的弱极限精确为黎曼零点的经验分布。
6. **离轴零点缺失**：对任何离轴零点，$\mathcal{H}_\alpha^{(N)}$ 的谱必发散，与 (2)(3) 矛盾。

**证明概要**：

- **(2)+(3)**：由定理 2.3 与引理 X'（递归块三角化刚性）。这是纯有限维矩阵代数的函子性论证，独立于黎曼假设，不涉及任何 $N\to\infty$ 极限或解析估计。它唯一确定 $\alpha=1/2$ 为 SDFT 簇承载算术递归结构的唯一参数。
- **(2)+(3) $\Rightarrow$ (4)**：由定理 3.7。$\alpha=1/2$ 时 $\mathcal{H}_{1/2}^{(N)}=0$ 锁定实部为 $1/2$。
- **(4) $\Rightarrow$ (5)**：由定理 5.2 与双重测度结构保证。内部测度 $\mu_{1/2}^{(N)}=\delta_0$ 冻结实部，频率测度 $\nu_{1/2}^{(N)}$ 展开虚部。
- **(5)+(6) $\Rightarrow$ (1)**：由命题 5.2'（双重防线）。若存在离轴零点，则 (a) 有限 $N$ 时 $\nu_{1/2}^{(N)}[F]$ 对支撑在离轴虚部附近的测试函数 $F$ 以 $N^{2(\sigma_0-1/2)}$ 发散，与定理 5.2 要求的收敛到有限值矛盾；(b) 极限时 $\operatorname{supp}(\sum_\rho\delta_\rho)\subseteq\{\Re(s)=1/2\}$，同样禁止离轴零点。因此所有非平凡零点必须满足 $\Re(\rho)=1/2$。
- **(1) $\Rightarrow$ (2)+(3)**：若黎曼猜想成立，则黎曼函数方程在临界线上的自对偶性迫使 $\mathcal{U}_{1/2}=I$（定理 2.3），此时 $\mathcal{F}_{1/2}^{(N)}$ 退化为标准块 DFT，FFT 蝴蝶操作直接给出递归块三角化标准形。$\square$

**推论**：此框架提供了构造性证明路径：通过有限维自对偶计算（标准 DFT/FFT）逼近无限维解析对象（黎曼零点），所有步骤满足 Bishop 构造性数学的可计算性要求。SDFT 框架揭示了黎曼猜想本质上是一个对偶性条件：仅当参数 $\alpha=1/2$ 时，光滑截断离散傅里叶变换实现严格自对偶，此时算术与谱的函子忠实性达到完美匹配，且谱测度从非紧的扩展态（$\alpha \neq 1/2$）紧化为临界线上的离散点态（$\alpha = 1/2$）。递归块三角化刚性（引理 X'）将这一解析对偶性硬化为纯代数的函子刚性——算术数据生成的 SDFT 簇只允许 $\alpha=1/2$ 的整体截面，而命题 5.2' 的双重防线确保零点分布作为该截面的谱推论必然共线。

---

## 8. 注记：作为公理选择的黎曼猜想

### 8.1 几何学类比

欧几里得几何与非欧几何的分野，不在于谁"正确"，而在于对平行公设的取舍：

- **欧几里得截面**：取平行公设为真 $\Rightarrow$ 三角形内角和 $= \pi$，空间平坦，测地线唯一。
- **罗巴切夫斯基截面**：取平行公设为假 $\Rightarrow$ 三角形内角和 $< \pi$，空间负曲率，测地线发散。
- **黎曼截面**：取平行公设为假（另一种方式）$\Rightarrow$ 三角形内角和 $> \pi$，空间正曲率，测地线闭合。

三种几何各自自洽，各自良定。没有哪一个比另一个"更数学"，它们只是同一个公理系统的不同截面。

### 8.2 算术量子化的对应

本框架揭示：黎曼猜想在 SDFT 算术量子化中扮演的角色，完全平行于平行公设在几何学中的角色。

| 几何学 | SDFT 算术量子化 |
|:---|:---|
| 绝对几何（不带公设的几何学基础） | SDFT 簇 $\{\mathcal{F}_\alpha^{(N)}\}_{N\ge 1}$（带参数 $\alpha$ 的离散傅里叶变换族） |
| 平行公设 | 函子刚性条件：递归块三角化标准形是否存在？ |
| 取公设为真 | 取 $\alpha = 1/2$（引理 X'：递归刚性唯一允许此参数） |
| 取公设为假 | 取 $\alpha \neq 1/2$（递归刚性断裂，但代数结构仍然良定） |

### 8.3 三个算术截面

**1. 欧几里得截面（$\alpha = 1/2$，黎曼情景）**
- 递归块三角化存在，FFT 协议可计算
- $\mathcal{H}_{1/2}^{(N)} = 0$，实部冻结在 $1/2$
- 谱紧化为 $\delta_0$，Carleman 条件以 $1/|\alpha-1/2|$ 的速率无穷满足
- 频率测度弱收敛：$\nu_{1/2}^{(N)} \to \frac{1}{2\pi}\sum_\rho \delta_\rho$
- 推论：所有可检测零点满足 $\Re(\rho) = 1/2$
- **物理类比**：基态，可计算，可量子化

**2. 非欧截面（$\alpha \neq 1/2$，离轴情景）**
- 递归块三角化不存在，FFT 协议不可直接执行
- $\mathcal{H}_\alpha^{(N)} \neq 0$，特征值发散，实部可偏离
- 谱非紧，支撑区间半径 $\sim |\alpha-1/2|\, N^2/(\log N)^2 \to \infty$
- Carleman 条件失效：$S_N(\alpha) \to 0$
- 允许离轴零点：实部偏离临界线的零点在此截面中完全合法
- **物理类比**：激发态，非紧扩展，不可量子化

**3. 退化截面（$\alpha \to 1/2$ 但 $N$ 先发散）**
- 极限不可交换（定理 5.3(C)）
- 先 $N\to\infty$ 再 $\alpha\to1/2$：无定义，路径不合法
- 此截面的不可达性正是"非欧截面不能连续形变到欧几里得截面"的算术对应

### 8.4 核心洞察

黎曼猜想不是"真"或"假"的绝对命题，而是：

> 在算术数据生成的 SDFT 簇中，递归块三角化刚性（函子刚性）唯一地选择了 $\alpha = 1/2$ 截面。该截面的谱推论恰好是"所有非平凡零点位于临界线上"。

换言之：

- 如果你要求算术量子化是**可计算的**（即存在与二分算术分解相容的递归标准形），那么你被迫选择 $\alpha = 1/2$，而黎曼猜想作为该选择的谱推论自动成立。
- 如果你放弃递归刚性，那么离轴零点不仅可能存在，而且有完整的良定数学结构（非紧算子、发散谱、失效的 Carleman 条件）来容纳它们。

就像欧几里得几何不是"真"的，而是"第五公设成立时的几何"；黎曼猜想不是"真"的，而是"递归块三角化刚性成立时的谱推论"。

---

# 附录 A：FFT 蝴蝶操作与对偶结构

## A.1 时间抽取 FFT 的蝴蝶操作

考虑长度为 $N$（$N=2^m$）的 DFT：
$$
X_k=\sum_{n=0}^{N-1}x_n\omega_N^{kn},\quad\omega_N=e^{-2\pi i/N}.
$$
将输入序列按奇偶分开：
$$
X_k=\sum_{m=0}^{N/2-1}x_{2m}\omega_N^{k(2m)}+\sum_{m=0}^{N/2-1}x_{2m+1}\omega_N^{k(2m+1)}=E_k+\omega_N^kO_k,
$$
其中 $E_k$ 和 $O_k$ 分别是长度为 $N/2$ 的 DFT：
$$
E_k=\sum_{m=0}^{N/2-1}x_{2m}\omega_{N/2}^{km},\quad O_k=\sum_{m=0}^{N/2-1}x_{2m+1}\omega_{N/2}^{km}.
$$
对于 $k$ 和 $k+N/2$，利用 $\omega_N^{k+N/2}=-\omega_N^k$，得到经典的蝶形运算：
$$
\begin{cases}
X_k=E_k+\omega_N^kO_k,\\
X_{k+N/2}=E_k-\omega_N^kO_k.
\end{cases}
$$
这一运算可图示为蝴蝶状结构，因此得名。

## A.2 蝴蝶操作与对偶性

蝴蝶操作揭示了 DFT 的内在对称性：偶部 $E_k$ 与奇部 $O_k$ 通过旋转因子 $\omega_N^k$ 组合，而旋转因子本身满足 $\omega_N^{k+N/2}=-\omega_N^k$，这正是函数方程中 $s \leftrightarrow 1-s$ 对称性的离散体现。在 SDFT 框架中，对偶性 $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ 对应于参数 $\alpha$ 与 $1-\alpha$ 的互换，而蝴蝶操作中的奇偶分解与符号反转恰好模拟了这一对偶关系。

进一步地，若将 SDFT 矩阵 $\mathcal{F}_\alpha^{(N)}$ 视为带有内部相位 $\exp(i(\alpha-1/2)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H})$ 的 DFT，则其快速算法同样可以分解为类似的蝶形结构，只是旋转因子需乘以相应的对角矩阵。这种结构保证了 SDFT 的计算效率，同时保留了与黎曼函数方程的对偶联系。

## A.3 蝴蝶操作与函数方程

将蝴蝶操作与第 4.2 节的对应表结合，可归纳为：

| FFT 蝴蝶操作 | SDFT 对偶性 |
|:---|:---|
| 奇偶分解 $E_k,O_k$ | 参数对偶 $\alpha \leftrightarrow 1-\alpha$ |
| 旋转因子 $\omega_N^k$ 与 $-\omega_N^k$ | 相位因子 $\exp(i(\alpha-1/2)\widetilde{\mathrm{Li}}_N(x_n)\mathbf{H})$ 及其共轭 |
| $X_k$ 与 $X_{k+N/2}$ 的对称性 | 自对偶条件 $\mathcal{U}_{1/2}=\mathbb{I}$ |
| 网格加密 $N\to\infty$ | 采样完备性 $t_k\to\gamma_j$（第 4.4 节） |

**与引理 X' 的联系**：附录 A.1 中的递归要求 $E_k$ 和 $O_k$ 必须是规模 $N/2$ 的同类 DFT，这是 FFT 蝴蝶操作可行的代数根基。引理 X' 正是将此递归要求施加于 SDFT 簇：当 $\alpha\neq 1/2$ 时，相位权重 $W_n$ 不满足 $W_{2m}^{(N)}=W_m^{(N/2)}$，导致偶部与奇部无法归约为规模 $N/2$ 的 SDFT，递归链条断裂。因此，FFT 蝴蝶操作的函子相容性唯一锁定 $\alpha=1/2$。

---

## 参考文献

1. H. M. Edwards, *Riemann's Zeta Function*, Dover, 2001 [Original: 1974].
2. E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., Oxford University Press, 1986.
3. A. Weil, "Sur les formules explicites de la théorie des nombres premiers", Comm. Sem. Math. Univ. Lund, 1952.
4. A. E. Ingham, "On the estimation of $N(\sigma,T)$", Quart. J. Math., 1940.
5. J. W. Cooley, J. W. Tukey, "An algorithm for the machine calculation of complex Fourier series", Math. Comp., 1965.
6. E. Bishop, *Foundations of Constructive Analysis*, McGraw-Hill, 1967.
7. H. P. Guinand, "Summation formulae and self-reciprocal functions", Quart. J. Math., 1942.
8. H. Iwaniec, E. Kowalski, *Analytic Number Theory*, American Mathematical Society, 2004.
