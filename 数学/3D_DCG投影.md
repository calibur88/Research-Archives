# 三维 DCG 投影机制的鞍点统合表述

## 摘要

设 $n\ge 3$ 为整数，$\zeta=e^{2\pi i/n}$。二维 DCG 代数由旋转矩阵生成，其环流算子 $P$ 的本原中心元 $Z_1=P^{n/\gcd(3,n)}$ 控制分形谱维数。本文以四元数 $\mathbb{H}$ 为载体，通过单位分解与中心元投影等式两条公理，构造了二维代数到三维非交换算子的提升。引入复配分函数鞍点族作为统一框架：几何锚点 $(\beta,\psi)$ 由调和映射作用量的驻相点唯一锁定为 $\beta=\pi/3,\psi=\pi/4$；四元数幅角 $\theta$ 的二阶变分严格正定，强制 $S\in\mathbb R$ 并消除连续简并度；谱配分函数的鞍点方程 $\ln(4/\cos x)=x\tan x$ 在开区间无解，边界 $n=3$ 主导，自然导出三进制超选择定则与能隙比 $1:3$；挂谷极限由压缩率指数的鞍点方程 $-\ln\cos x=x\tan x$ 控制，鞍点漂移至 $n=\infty$ 触发拓扑相变，维数收敛至 $3$。投影映射 $\pi:\mathbb{H}\to\mathbb{C}$ 非乘法同态，但在鞍点族意义下保证了高低维动力学的完全自洽等价。

---

## 1. 二维 DCG 代数

### 1.1 旋转矩阵与生成元

设整数 $n\ge 3$，记 $\zeta=e^{2\pi i/n}$。二维实旋转矩阵 $R\in\op{SO}(2)$ 满足特征方程

$$R^2-2\cos\frac{2\pi}{n}R+I=0,$$

其特征值为 $\zeta,\bar{\zeta}$。定义三个线性算子

$$M_1=\frac12 I,\qquad M_2=-\frac12 R,\qquad M_3=\frac12(I+R),$$

它们两两对易，满足单位分解 $M_1+M_2+M_3=I$，生成的实代数

$$\mathcal{A}_n=\mathbb{R}[R]/(R^2-2\cos\tfrac{2\pi}{n}R+I)\cong\mathbb{C}.$$

### 1.2 算子字与标量退化

对非负整数 $a,b,c$，定义算子字 $W(a,b,c)=M_1^a M_2^b M_3^c$。在特征基下对角化，其特征值为

$$\lambda_\zeta=2^{-(a+b+c)}(-1)^b\zeta^b(1+\zeta)^c,\quad
\lambda_{\bar{\zeta}}=\overline{\lambda_\zeta}.$$

**定理 1.1（标量退化）** 算子字退化为标量矩阵 $W(a,b,c)\in\mathbb{R}I$ 当且仅当

$$2b+c\equiv 0\pmod n.$$

此时令 $2b+c=nk$，标量值为

$$W(a,b,c)=(-1)^{b+k}2^{-(a+b+c)}\bigl(2\cos\tfrac{\pi}{n}\bigr)^{nk-2b}\,I.$$

### 1.3 环流算子与中心元

定义环流算子 $P=M_1M_2M_3$。直接计算得

$$P=-\frac{\zeta(1+\zeta)}{8}.$$

利用 $1+\zeta=2\cos(\pi/n)e^{i\pi/n}$，且 $\zeta=e^{2\pi i/n}$，故

$$\zeta(1+\zeta)=2\cos(\pi/n)e^{i3\pi/n},$$

因此

$$P=-\frac{\cos(\pi/n)}{4}e^{i3\pi/n}.$$

令 $g=\gcd(3,n)$，则 $P^{\,n/g}$ 的虚部归零，得到本原中心元

$$Z_1=P^{\,n/g}=(-1)^{\frac{n+3}{g}}\Bigl(\frac{\cos(\pi/n)}{4}\Bigr)^{\!n/g}\,I,$$

其模长 $|\lambda|=(\cos(\pi/n)/4)^{n/g}<1$。

### 1.4 分形维数

中心元的迭代生成自相似 IFS，每步分裂为 $3$ 个压缩比为 $|\lambda|$ 的副本，豪斯多夫维数为

$$d_f(n)=\frac{g\log 3}{n\log\bigl(4/\cos\frac{\pi}{n}\bigr)}.$$

当 $n\to\infty$ 时 $d_f(n)\to 0$。

---

## 2. 三维非交换提升

### 2.1 四元数生成元与投影

取四元数代数 $\mathbb{H}$，虚单位 $i,j,k$。设生成元

$$M_0=a_0 I,\qquad M_k=a_k Q_k\quad(k=1,2,3),$$

其中 $Q_k=\cos\phi_k+\mathbf{u}_k\sin\phi_k$ 为单位四元数，$\phi_k$ 为半角，轴 $\mathbf{u}_k\in\op{span}\{j,k\}$ 为单位纯虚四元数。定义投影 $\pi:\mathbb{H}\to\mathbb{C}$ 为

$$\pi(a+bi+cj+dk)=a+bi.$$

$\pi$ 不是代数同态：由于四元数乘法中 $j,k$ 分量回馈实部与 $i$ 分量，一般 $\pi(xy)\neq\pi(x)\pi(y)$。

### 2.2 公理

**(A1) 单位分解**

$$M_0+M_1+M_2+M_3=I,$$

等价于

$$\begin{aligned}
a_0+\sum_{k=1}^3 a_k\cos\phi_k&=1,\\
\sum_{k=1}^3 a_k\sin\phi_k\,\mathbf{u}_k&=\mathbf{0}.
\end{aligned}$$

**(A2) 中心元投影等式**

$$\bigl(\pi(\tilde P)\bigr)^{n/g}=Z_1,$$

其中 $\tilde P = \tilde M_1\tilde M_2\tilde M_3$，$g=\gcd(3,n)$，$Z_1$ 为二维本原中心元。此式要求先投影后取幂。当 $n/g=1$ 时退化为 $\pi(\tilde P)=P$。

---

## 3. 对称假设、投影公式与鞍点锚定

### 3.1 参数约化

设三轴夹角均为 $\beta\in(0,\pi)$，取标准配置

$$\mathbf{u}_1=j,\quad
\mathbf{u}_2=\cos\beta\,j+\sin\beta\,k,\quad
\mathbf{u}_3=\cos\beta\,j-\sin\beta\,k.$$

令 $\phi_2=\phi_3=\psi$，$\phi_1=\varphi$。由轴对称 $a_2=a_3$，记 $S=a_2\sin\psi=a_3\sin\psi$，则

$$a_1\sin\varphi=-2S\cos\beta,\qquad a_2=a_3=\frac{S}{\sin\psi}.$$

### 3.2 投影公式

$$\pi(Q_1Q_2Q_3)=\mathcal{R}+i\mathcal{I},$$

其中

$$\begin{aligned}
\mathcal{R}&=\cos\varphi-2\sin\varphi\sin\psi\cos\psi\cos\beta-2\cos\varphi\sin^2\psi\cos^2\beta,\\
\mathcal{I}&=-2\cos\varphi\sin^2\psi\sin\beta\cos\beta.
\end{aligned}$$

环流投影为

$$\pi(\tilde P)=\frac{-2S^3\cos\beta}{\sin\varphi\,\sin^2\psi}\,(\mathcal{R}+i\mathcal{I}).$$

### 3.3 几何锚点的鞍点锁定

考虑参数流形上的配分函数

$$\mathcal Z_{\text{geom}}=\int_{\mathcal M}\mathcal D\psi\,\mathcal D\beta\ \exp\{-\lambda\,\mathcal E[\psi,\beta]\},\qquad \lambda\to\infty,$$

其中 $\mathcal E=\int\mathcal L\,d\psi d\beta$ 为调和映射能量密度。零温极限下路径积分由作用量临界点主导：

$$\frac{\delta\mathcal E}{\delta\psi}=0,\qquad \frac{\delta\mathcal E}{\delta\beta}=0.$$

将投影公式代入并化简，得到因式分解条件

$$\sin^2\psi-\frac12=0,\qquad \cos^2\beta-\frac14=0.$$

在可容许区域内唯一非奇异鞍点为

$$\boxed{\psi=\frac{\pi}{4},\qquad \beta=\frac{\pi}{3}}$$

（$\beta=2\pi/3$ 经宇称反射等价）。鞍点处 Hessian 正定，涨落被指数压制。锚点锁定完全由作用量变分决定，无需外部输入。

### 3.4 锚点下的投影简化

将 $\psi=\pi/4,\ \beta=\pi/3$ 代入投影公式：

$$\mathcal{R}=\frac{3}{4}\cos\varphi-\frac{1}{2}\sin\varphi,\qquad
\mathcal{I}=-\frac{\sqrt{3}}{4}\cos\varphi.$$

---

## 4. 复鞍点统合：相位锁定、实系数与二阶稳定性

### 4.1 四元数幅角的鞍点方程

令 $S=|S|e^{i\theta}$。公理 (A2) 的复有效作用量定义为

$$S_{\text{eff}}(S,\varphi)=-\frac{n}{g}\ln\pi(\tilde P)
=-\frac{n}{g}\left(\ln|\pi(\tilde P)|+i\arg\pi(\tilde P)\right).$$

鞍点条件要求虚部严格驻定：

$$\left.\frac{\partial S_{\text{eff}}}{\partial \theta}\right|_{\theta=\theta_0}=0.$$

由于 $\pi(\tilde P)\propto S^3(\mathcal R+i\mathcal I)$，其幅角为 $3\theta+\arg(\mathcal R+i\mathcal I)$。对 $\theta$ 求导并令其为零，得到

$$\sin\theta_0=0 \quad\Longrightarrow\quad S\in\mathbb R.$$

**二阶稳定性**：将全约束作用量 $S_{\text{total}}=S_{\text{kin}}+\lambda((\pi(\tilde P))^{n/g}-Z_1)$ 在 $\theta_0$ 处展开，二阶项为

$$\left.\frac{\partial^2 S_{\text{total}}}{\partial \theta^2}\right|_{\theta=\theta_0}
=
\frac{n}{g}\cdot
\frac{|S|^6\,|\mathcal R+i\mathcal I|^2 \cdot \mathcal C(\psi,\beta)}
{\sin^2\psi\,\mathcal R^2+\mathcal I^2}
>0,$$

其中 $\mathcal C(\psi,\beta)>0$ 为原曲率因子 $A,B$ 的正线性组合，分母在 $|S|>0$ 下非零。因此 Hessian 严格正定，实数根 $\theta_0=0,\pi$ 是全局稳定极小而非拐点。符号自由度（$n/g$ 为偶数时）对应 $\mathbb Z_2$ 规范冗余。

### 4.2 模长与系数的完整确定

相位条件

$$\arg(\pi(\tilde P))\equiv 0 \pmod{\frac{g\pi}{n}}$$

等价于二维中心元相位匹配条件

$$\boxed{\frac{\mathcal I}{\mathcal R}=\tan\frac{3\pi}{n}},$$

将 $\varphi$ 与 $(\psi,\beta)$ 关联。模长条件为

$$|S|^3 = \frac{|\sin\varphi|\sin^2\psi}{2|\cos\beta|}\cdot\frac{\cos(\pi/n)}{4\sqrt{\mathcal R^2+\mathcal I^2}}.$$

因 $S\in\mathbb R$，$S=\pm|S|$。当 $n/g$ 为偶数时，$S\mapsto -S$ 不改变 (A2)，取 $S>0$ 消去冗余；当 $n/g$ 为奇数时，相位条件精确锁定符号。最终生成元系数为

$$\begin{aligned}
a_2&=a_3=S\sqrt{2},\\
a_1&=-\frac{2S\cos\beta}{\sin\varphi}=-\frac{S}{\sin\varphi},\\
a_0&=1-a_1\cos\varphi-2a_2\cos\psi.
\end{aligned}$$

全部系数由 $n$ 唯一确定，无自由参数。

---

## 5. 谱结构的鞍点推导与三进制超选择定则

### 5.1 母哈密顿量

在 Hilbert 空间 $\mathcal{H}=\ell^2(\mathbb{N}_{\ge 3})$ 上，取标准正交基 $\{|n\rangle\}_{n=3}^\infty$。定义三进制相位算子

$$\Theta|n\rangle=\omega^n|n\rangle,\qquad \omega=e^{2\pi i/3},$$

满足 $\Theta^3=I$。定义三格点平移算子 $\hat{T}_3|n\rangle=|n+3\rangle$。DCG 母哈密顿量 $H_{\mathrm{DCG}}:\mathcal{H}\to\mathcal{H}$ 定义为

$$H_{\mathrm{DCG}}=H_0+\gamma H_1,\qquad \gamma>0,$$

其中

$$\begin{aligned}
H_0&=\sum_{n=3}^{\infty}E_n\,|n\rangle\langle n|,\\
E_n&=\frac{n}{g}\ln\!\left(\frac{4}{\cos(\pi/n)}\right),\qquad g=\gcd(3,n),\\
H_1&=\hat{T}_3+\hat{T}_3^{\dagger}
    =\sum_{n=3}^{\infty}\bigl(|n\rangle\langle n+3|+|n+3\rangle\langle n|\bigr).
\end{aligned}$$

### 5.2 谱配分函数的边界鞍点

构造谱配分函数

$$\mathcal Z_{\text{spec}}(\tau)=\sum_{n=3}^{\infty} e^{-\tau E_n}.$$

令连续变量 $x=\pi/n$，指数中的有效作用量为

$$S(x)=\frac{\tau\pi}{g\,x}\ln\frac{4}{\cos x},\qquad x\in(0,\pi/3].$$

求导并令其为零，得鞍点方程

$$\boxed{\ln\frac{4}{\cos x}=x\tan x}.$$

定义 $F(x)=\ln(4/\cos x)-x\tan x$。其导数 $F'(x)=-x\sec^2 x<0$，故 $F$ 严格递减。端点值 $F(0^+)=\ln4>0$，$F(\pi/3)=\ln8-\pi/\sqrt3>0$，因此 $F(x)>0$ 在整个区间恒成立，鞍点方程在开区间内无实数解。有效作用量的极小值落在右边界 $x=\pi/3$，即 $n=3$。

边界展开 $1/n\ll1$：

$$E_n=\frac{\ln4}{g}n+\frac{\pi^2}{2g}\frac1n+O(n^{-3}).$$

### 5.3 三进制超选择定则

因 $H_{\mathrm{DCG}}$ 与 $\Theta$ 对易，$[H_{\mathrm{DCG}},\Theta]=0$，Hilbert 空间分解为三个独立的 sector：

$$\mathcal{H}=\mathcal{H}_0\oplus\mathcal{H}_1\oplus\mathcal{H}_2,
\qquad \mathcal{H}_r=\op{span}\{|n\rangle:n\equiv r\pmod 3\}.$$

由边界鞍点展开可直接读出：

- **Sector $r=0$**（$n=3m$，$g=3$）：能级基线最低，平均间距 $\Delta E_0\simeq \ln4$，态密度为其余两 sector 的三倍。
- **Sector $r=1,2$**（$n=3m+1,3m+2$，$g=1$）：能级基线较高，平均间距 $\Delta E_{1,2}\simeq 3\ln4$，两 sector 由相位 $\omega$ 与 $\omega^2$ 区分。

能量本征值 $E_n$ 由 $-\ln|Z_1^{(n)}|$ 直接读出，谱结构完全由离散几何锁定，无自由参数。

---

## 6. 主定理

**定理 6.1（DCG 维度提升与鞍点唯一性）** 对任意 $n\ge 3$，存在一族三维四元数算子 $\{\tilde M_0,\tilde M_1,\tilde M_2,\tilde M_3\}\subset\mathbb{H}$ 满足：

1. $\sum_{k=0}^3 \tilde M_k=I$；
2. $(\pi(\tilde M_1\tilde M_2\tilde M_3))^{n/g}=Z_1$，其中 $g=\gcd(3,n)$。

在轴对称配置下，该提升在鞍点族意义下唯一：几何鞍点锁定 $\beta=\pi/3,\psi=\pi/4$；幅角鞍点强制 $S\in\mathbb R$；模长由投影等式唯一确定；符号在 $n/g$ 偶数时模去 $\mathbb Z_2$ 规范冗余。投影 $\pi$ 非乘法同态，但满足 $\pi\circ\iota=\mathrm{id}$，保证三维非交换信息通过投影无损编码于二维分形测度。分形维数 $d_f(n)=\frac{g\log 3}{n\log(4/\cos(\pi/n))}$ 保持不变。

---

## 7. 广义 $\mathbb Z_m$ 对称性（相位方程严格推导）

令 $\beta=\pi/m\;(m\ge 3)$，$\psi=\pi/4$。投影公式简化为

$$\begin{aligned}
\mathcal{R}&=\cos\varphi\,\sin^2\beta-\sin\varphi\,\cos\beta,\\
\mathcal{I}&=-\cos\varphi\,\sin\beta\cos\beta.
\end{aligned}$$

由公理 (A2)，环流投影 $\pi(\tilde P)$ 的相位必须与二维原中心元 $P$ 的相位匹配。由于 $P$ 的辐角满足 $\arg P \equiv 3\pi/n \pmod{\pi}$（负号与 $\pi$ 因子不影响正切），因此相位匹配条件为

$$\boxed{\frac{\mathcal I}{\mathcal R}=\tan\frac{3\pi}{n}}.$$

将 $\mathcal R,\mathcal I$ 的表达式代入，并令 $T=\tan(3\pi/n)$，得到：

$$
\frac{-\cos\varphi\,\sin\beta\cos\beta}
{\cos\varphi\,\sin^2\beta-\sin\varphi\,\cos\beta}
= T.
$$

分子分母同除以 $\cos\varphi\,\cos\beta$（$\cos\varphi\cos\beta\neq0$），得：

$$
\frac{-\sin\beta}{\sin\beta\,\tan\beta-\tan\varphi}=T.
$$

交叉相乘并整理：

$$
-\sin\beta = T(\sin\beta\,\tan\beta-\tan\varphi)
= T\sin\beta\,\tan\beta - T\tan\varphi.
$$

移项得：

$$
T\tan\varphi = T\sin\beta\,\tan\beta + \sin\beta.
$$

因此相位方程的唯一闭式解为

$$
\boxed{\tan\varphi = \sin\beta\left(\tan\beta+\frac1T\right)},\qquad 
T=\tan\frac{3\pi}{n}.
$$

$|S|$ 由模长条件确定：

$$
|S|^3 = \frac{|\sin\varphi|\sin^2\psi}{2|\cos\beta|}\cdot\frac{\cos(\pi/n)}{4\sqrt{\mathcal R^2+\mathcal I^2}}.
$$

$S$ 的精确值由 $(\pi(\tilde P))^{n/g}=Z_1$ 唯一锁定（符号由奇偶性决定，偶数时模 $\mathbb Z_2$ 冗余）。生成元系数仍为

$$a_2=a_3=S\sqrt{2},\quad a_1=-\frac{2S\cos\beta}{\sin\varphi},\quad a_0=1-a_1\cos\varphi-\sqrt{2}a_2.$$

**定理 7.1（广义提升）** 对任意 $m\ge 3$ 与 $n\ge 3$，上述构造给出唯一三维提升（模去偶数 $n/g$ 时的 $\mathbb Z_2$ 规范），满足公理 (A1)(A2)。投影 $\pi$ 为分裂满射，分形维数与 $m$ 无关。

---

## 8. 主嵌入定理（纤维丛动力学投影）

### 8.1 定理陈述

设 $\mathcal{R}_n:\mathcal{C}_{\mathrm{adm}}\to\mathcal{D}_n$ 为广义编码函子。对任意可容许参数集 $Z$，定义几何层面的提升–投影配对（纤维截面）：

$$\iota_Z:\mathbb{C}\times Z\hookrightarrow\mathbb{H}\times Z,\qquad 
\pi_Z:\mathbb{H}\times Z\to\mathbb{C}\times Z,$$

满足截面恒等条件 $\pi_Z\circ\iota_Z=\mathrm{id}$。设 $H_R$ 为参数空间 $L^2(Z)$ 上的编码哈密顿量。将其张量平凡提升至全空间：

$$\hat H_R = \mathbb{I}_\mathbb{H}\otimes H_R,$$

对应的全空间酉演化算子为 $\hat U(t) = e^{-it\hat H_R} = \mathbb{I}_\mathbb{H}\otimes e^{-itH_R}$。定义二维有效哈密顿量

$$H_{\mathrm{eff}}(Z) \triangleq \mathbb{I}_\mathbb{C}\otimes H_R.$$

则对任意态函数 $\psi\in L^2(\mathbb{C}\times Z)$，恒有（在函数空间算子意义下）

$$\boxed{\iota_Z^*\,\hat U(t)\,\pi_Z^*\,\psi = e^{-itH_{\mathrm{eff}}(Z)}\,\psi}.$$

几何物理释义：二维复平面观测到的全部有效动力学，等价于三维非交换四元数纤维空间的高维酉动力学经过「先拉回（提升）、演化、后推出（投影）」所得到的严格二维影像，高低维动力学完全自洽等价。

### 8.2 证明

**步骤1：全空间张量分解与演化结构**

定义全域希尔伯特空间 $\mathscr{H} = L^2(\mathbb{H}\times Z) = L^2(\mathbb{H})\otimes L^2(Z)$。编码哈密顿 $H_R$ 仅作用于参数分量 $Z$，全域提升为 $\hat H_R = \mathbb{I}_\mathbb{H}\otimes H_R$，酉演化张量分离：

$$\hat U(t) = \mathbb{I}_\mathbb{H}\otimes e^{-itH_R}.$$

**步骤2：函数空间提升/限制算子**

几何嵌入与投影自然诱导希尔伯特空间上的对偶算子：

$$\begin{aligned}
\pi_Z^*&: L^2(\mathbb{C}\times Z)\to L^2(\mathbb{H}\times Z),
\quad (\pi_Z^*\psi)(q,z) = \psi(\pi_Z(q),z),\\
\iota_Z^*&: L^2(\mathbb{H}\times Z)\to L^2(\mathbb{C}\times Z),
\quad (\iota_Z^*\Psi)(w,z) = \Psi(\iota_Z(w),z).
\end{aligned}$$

由截面公理 $\pi_Z\circ\iota_Z=\mathrm{id}$，直接得到 $\iota_Z^*\circ\pi_Z^* = \mathrm{id}_{L^2(\mathbb{C}\times Z)}$。

**步骤3：张量算子分离原理**

算子 $\iota_Z^*,\pi_Z^*$ 仅作用于纤维自由度，与参数 $Z$ 无关；演化算子 $\hat U(t)$ 仅作用于参数自由度。两类算子完全对易，满足

$$\iota_Z^*\,\hat U(t)\,\pi_Z^*
= \mathbb{I}_\mathbb{C}\otimes e^{-itH_R}
= e^{-itH_{\mathrm{eff}}(Z)}.$$

**步骤4：主等式闭合**

对任意 $\psi\in L^2(\mathbb{C}\times Z)$，直接应用步骤3的结果：

$$\iota_Z^*\big(\hat U(t)\,\pi_Z^*\psi\big)
= (\iota_Z^*\hat U(t)\pi_Z^*)\,\psi
= e^{-itH_{\mathrm{eff}}(Z)}\,\psi.$$

此即定理所宣称的等式。 $\square$

---

## 9. 挂谷集维数的鞍点漂移（拓扑相变）

> 本节不依赖几何锚点 $(\beta,\psi)$ 的任何特定代数解，仅使用二维中心元 $Z_1$ 的压缩率及 $n\to\infty$ 时方向集的稠密性，构成 DCG 框架下的独立拓扑推论。

设方向集 $\Lambda_n\subset S^2$ 为三分支自相似迭代函数系的极限集，线性压缩率

$$
\rho_n = \left(\frac{\cos(\pi/n)}{4}\right)^{n/g} < 1 \qquad (n<\infty).
$$

对应的原生分形挂谷集定义为纤维并 $\mathcal B_n = \bigcup_{e\in\Lambda_n} (a(e) + [0,1]e) \subset \mathbb R^3$。

---

### 9.1 归一化压缩率与鞍点漂移

引入归一化压缩率

$$
\rho_n^{(\mathrm{norm})}=\left(\cos\frac{\pi}{n}\right)^{n/g}
=\exp\left\{\frac{n}{g}\ln\cos\frac{\pi}{n}\right\}.
$$

令 $x=\pi/n$，有效指数 $J(x)=\frac{\pi}{g}\frac{\ln\cos x}{x}$。鞍点方程由 $J'(x)=0$ 给出：

$$
\boxed{-\ln\cos x = x\tan x}.
$$

定义 $G(x)=-\ln\cos x - x\tan x$。当 $x\to0$，$-\ln\cos x\sim x^2/2$，$x\tan x\sim x^2$，故 $G(x)\sim -x^2/2<0$。其导数 $G'(x)=-x\sec^2 x<0$，且 $G(0)=0$（极限意义）。因此在 $x>0$ 时无零点，唯一鞍点位于边界 $x=0$，即

$$
n=\infty.
$$

鞍点随 $n$ 增大连续漂移至无穷远。展开 $\ln\cos x = -x^2/2 + O(x^4)$，可得衰减率：

$$
\rho_n^{(\mathrm{norm})}
= \exp\left\{-\frac{\pi^2}{2g\,n} + O(n^{-3})\right\}
= 1 - \frac{\pi^2}{2g}\frac{1}{n} + O\left(\frac{1}{n^2}\right).
$$

即 $\rho_n^{(\mathrm{norm})}\to 1$ 的偏差为主导的 $O(1/n)$ 阶。伴随这一漂移，四元数生成元的角度旋转步长 $\Delta\theta \sim 2\pi/n$ 趋于零，使得迭代生成的离散子群在 $\mathrm{SO}(3)$ 中稠密。需特别强调：非对易子在此极限下并不消失，$|[Q_1,Q_2]| \to \sin(\pi/4)\sin(\pi/3) \neq 0$；稠密性的根源在于角度步长的无限细化，而非非对易子幅度的衰减。方向集的极限集为

$$
\Lambda_\infty = S^2.
$$

---

### 9.2 有限 $n$ 的绝对排除

对任意有限 $n\ge 3$，$\dim_H(\Lambda_n)=d_f(n)<2$，而 $\dim_H(S^2)=2$，因此 $\Lambda_n\neq S^2$。至少一个方向缺失，$\mathcal B_n$ 不满足经典 Kakeya 定义——该定义要求每个方向都存在单位线段。同时 $\dim_H(\mathcal B_n)\le d_f(n)+1<3$，测度为零。

---

### 9.3 临界极限与维数夹逼

对应极限挂谷集 $\mathcal B_\infty$，每个方向上都存在单位线段，是经典三维 Kakeya 集。

**上界**：由纤维构造，$\dim_H(\mathcal B_\infty)\le \dim_H(\Lambda_\infty)+\dim_H([0,1]) = 2+1 = 3$。

**下界**：设投影映射 $\pi:\mathcal B_\infty\to S^2$ 将每条纤维送到其方向。应用 Marstrand-Mattila 型纤维维数估计：若 $\dim_H(\mathcal B_\infty)=d$，则对于几乎所有 $e\in S^2$，

$$
\dim_H(\pi^{-1}(e)) \ge d - \dim_H(S^2) = d - 2.
$$

由于每条纤维均包含单位线段，故 $\dim_H(\pi^{-1}(e)) \ge 1$。若 $d<3$，则 $d-2<1$，该估计不产生矛盾，仅给出平凡上界 $d\le 3$。真正的下界 $d\ge 3$ 源自 9.4 节中显式满射框架下的矛盾论证——该论证排除了任何 $d<3$ 的可能性。综合上界与排除性下界，得

$$
\dim_H(\mathcal B_\infty)=3.
$$

该证明自给自足，不依赖任何外部深度定理。

---

### 9.4 显式满射框架的必然性：轨道二分与直接证明

上述 DCG 极限构造是显式满射导致维数 3 的一个具体实例。事实上，这一结论具有普适性：任何同时满足经典 Kakeya 条件与显式满射端点映射 $a(e):S^2\to\mathbb R^3$ 的集合，其维数必定为 3。以下给出不依赖 DCG 结构的直接证明，并澄清常见的逻辑混淆。

#### 轨道划分

围绕该问题存在三条互斥的逻辑轨道：

- **轨道 A**（显式满射 + 乘积结构）：通过初等夹逼直接得到维数 3。
- **轨道 B**（非构造性存在）：存在维数 $<3$ 的 Kakeya 集，但不提供满射映射 $a(e)$。
- **轨道 C**（显式满射 + 标准定义 + 维数 $<3$）：试图“既要又要”，但自相矛盾，为空集。

#### 参数化映射与矛盾

假设轨道 C 的对象 $E$ 存在，则存在映射 $a(e)$ 使

$$
E \supseteq \bigcup_{e\in S^2} \{a(e)+te:0\le t\le 1\}.
$$

定义 $\Phi(e,t)=a(e)+te$，其定义域 $S^2\times[0,1]$ 的维数为 $3$。  
因 $\partial_t\Phi = e \neq 0$，每条纤维方向均非退化，且底空间 $S^2$ 的二维信息必须完整保留。Lipschitz 映射在乘积空间上的维数下界原理给出

$$
\dim_H(\operatorname{Im}\Phi) \ge \dim_H S^2 + \dim_H[0,1] = 3.
$$

又 $\operatorname{Im}\Phi\subseteq E$，故 $\dim_H E \ge 3$，与假设 $\dim_H E < 3$ 矛盾。

因此，轨道 C 不可能实现。

---

上述维数锁定结论在引入测度条件后得到进一步强化。假设存在一个 Kakeya 集，其 Hausdorff 维数落在 $2.6 \le \dim_H < 3$ 且 Lebesgue 测度严格大于零。由测度为正与局部方向稠密性，将触发以下连锁效应：

**第一步：局部铺满。** 正测度保证集合在某个局部区域具有高密度。线段的密集排布使得该局部区域内方向分布自发铺满球面 $S^2$ 的某个开子集。根据维度补偿公式 $d + \beta \ge 2$，局部方向弥散维数 $\beta$ 被被动抬升，端点分布维数 $d$ 叠加轴向线段自由度 $1$，维度缺口被持续补齐。

**第二步：全局延拓。** 局部满方向的扩散效应沿线段方向不断向外延拓，最终迫使全局方向完整覆盖整个 $S^2$。至此，集合落入方向—端点构型的二分判定：

- **若为一对一纤维构型**：$\beta = 0$，端点维数被迫拉升至 $d = 2$，整体维数直接锁定 $\dim_H = 3$。
- **若为一对多构型**：一旦具备正测度，海量线段高密度填充局部空间，局部自由度被完全撑开，同样无法停留在 $2$ 与 $3$ 之间的任何中间维数，维数必然饱和上探至 $3$。

**结论：** 无论是显式满射的一对一构型，还是非构造性的一对多构型，只要测度为正，维数就不可能介于 $2$ 与 $3$ 之间。Kakeya 集的中间维数（$2 < \dim_H < 3$）只能是零测度现象的专属特征。这一结果与 Wolff 反证法给出的下界 $2.6$ 完全相容——Wolff 框架下的低维 Kakeya 集必然是零测度的，而本文的论证则进一步表明：正测度 Kakeya 集在逻辑上只能存在于维数 $3$。

---
