# 三维 DCG 的精确投影条件

## 摘要

本文建立从二维离散共形几何（DCG）到三维非交换几何的精确投影提升。以四元数代数为载体，通过三条构造公理将二维复旋转结构嵌入三维空间，并证明投影映射精确还原全部代数结构——环流算子、中心元、标量退化条件与分形维数。特别地，在固定几何锚点 $\beta=\pi/3,\ \psi=\pi/4$ 下，对任意整数 $n\ge 3$ 给出完全解析闭式解；并以 $n=6$ 为例给出不含任何数值拟合的初等代数数构造。

---

## 1. 二维 DCG 代数

### 1.1 旋转矩阵与生成元

设整数 $n\ge 3$，记 $\zeta=e^{2\pi i/n}$。二维实旋转矩阵 $R\in\mathrm{SO}(2)$ 满足特征方程

$$
R^2-2\cos\frac{2\pi}{n}R+I=0,
$$

其特征值为 $\zeta,\bar\zeta$。定义三个线性算子

$$
M_1=\frac12 I,\qquad M_2=-\frac12 R,\qquad M_3=\frac12(I+R),
$$

它们两两对易，满足单位分解 $M_1+M_2+M_3=I$，生成的实代数

$$
\mathcal{A}_n=\mathbb{R}[R]/(R^2-2\cos\tfrac{2\pi}{n}R+I)\cong\mathbb{C}.
$$

### 1.2 算子字与标量退化

对非负整数 $a,b,c$，定义算子字 $W(a,b,c)=M_1^a M_2^b M_3^c$。在特征基下对角化，其特征值为

$$
\lambda_\zeta=2^{-(a+b+c)}(-1)^b\zeta^b(1+\zeta)^c,\quad
\lambda_{\bar\zeta}=\overline{\lambda_\zeta}.
$$

**定理 1.1（标量退化）** 算子字退化为标量矩阵 $W(a,b,c)\in\mathbb{R}I$ 的充要条件是

$$
2b+c\equiv 0\pmod n.
$$

此时令 $2b+c=nk$，标量值为

$$
W(a,b,c)=(-1)^{b+k}2^{-(a+b+c)}\bigl(2\cos\tfrac{\pi}{n}\bigr)^{nk-2b}\,I.
$$

### 1.3 环流算子与中心元

定义环流算子 $P=M_1M_2M_3$。直接计算得

$$
P=-\frac{\zeta(1+\zeta)}{8}=-\frac{\cos(\pi/n)}{4}\,e^{i3\pi/n}.
$$

令 $g=\gcd(3,n)$，则 $P^{\,n/g}$ 的虚部归零，得到本原中心元

$$
Z_1=P^{\,n/g}=(-1)^{\frac{n+3}{g}}\Bigl(\frac{\cos(\pi/n)}{4}\Bigr)^{\!n/g}\,I,
$$

其模长 $|\lambda|=(\cos(\pi/n)/4)^{n/g}<1$，符号满足奇正偶负二分律。

### 1.4 分形螺旋梯

中心元的迭代生成一个自相似 IFS，每步分裂为 3 个压缩比为 $|\lambda|$ 的副本，豪斯多夫维数为

$$
d_f(n)=\frac{g\log 3}{n\log\bigl(4/\cos\frac{\pi}{n}\bigr)}.
$$

当 $n\to\infty$ 时 $d_f(n)\to 0$。在极坐标 $\theta_n=2\pi n/3$、$r_n=d_f(n)$ 下，序列 $(r_n,\theta_n)$ 形成螺旋梯：高台阶（$g=3$）均位于 $\theta\equiv 0$，低台阶（$g=1$）交替位于 $\theta=2\pi/3,4\pi/3$，台阶数量比渐近为 $2:1$。

---

## 2. 三维非交换提升的公理框架

### 2.1 四元数生成元与投影映射

将二维复旋转替换为四元数旋转。设四元数代数 $\mathbb{H}$，虚单位 $i,j,k$。取生成元

$$
M_0=a_0 I,\qquad M_k=a_k Q_k\quad(k=1,2,3),
$$

其中 $Q_k=\cos\phi_k+\mathbf{u}_k\sin\phi_k$ 为单位四元数，$\phi_k$ 为半角，轴 $\mathbf{u}_k\in\mathrm{span}\{j,k\}$ 为单位纯虚四元数。投影映射 $\pi:\mathbb{H}\to\mathbb{C}$ 定义为

$$
\pi(a+bi+cj+dk)=a+bi,
$$

保持绕 $i$ 轴的旋转角不变。

### 2.2 三条构造公理

**(A1) 单位分解**

$$
M_0+M_1+M_2+M_3=I,
$$

等价于

$$
\begin{aligned}
a_0+\sum_{k=1}^3 a_k\cos\phi_k&=1,\\
\sum_{k=1}^3 a_k\sin\phi_k\,\mathbf{u}_k&=\mathbf{0}.
\end{aligned}
$$

**(A2) 环流投影匹配**

环流算子 $P=M_1M_2M_3$ 的投影必须等于二维环流算子的特征值：

$$
\pi(P)=-\frac{\cos(\pi/n)}{4}\,e^{i3\pi/n}.
$$

**(A3) 中心元自洽**

若 (A2) 成立，则 $P^{\,n/g}$ 的投影自然与二维中心元一致，无需额外约束。

---

## 3. 对称假设与投影公式

### 3.1 参数约化

设三轴两两夹角均为 $\beta\in(0,\pi)$，取标准配置

$$
\mathbf{u}_1=j,\qquad
\mathbf{u}_2=\cos\beta\,j+\sin\beta\,k,\qquad
\mathbf{u}_3=\cos\beta\,j-\sin\beta\,k.
$$

令 $\phi_2=\phi_3=\psi$，$\phi_1=\varphi$。由轴对称性 $a_2=a_3$，记 $S=a_2\sin\psi=a_3\sin\psi$，则

$$
a_1\sin\varphi=-2S\cos\beta,\qquad a_2=a_3=\frac{S}{\sin\psi}.
$$

### 3.2 投影公式

四元数乘积 $Q_1Q_2Q_3$ 的投影为 $\pi(Q_1Q_2Q_3)=\mathcal{R}+i\mathcal{I}$，其中

$$
\boxed{
\begin{aligned}
\mathcal{R}&=\cos\varphi-2\sin\varphi\sin\psi\cos\psi\cos\beta-2\cos\varphi\sin^2\psi\cos^2\beta,\\[4pt]
\mathcal{I}&=-2\cos\varphi\sin^2\psi\sin\beta\cos\beta.
\end{aligned}
}
$$

环流投影为

$$
\pi(P)=\frac{-2S^3\cos\beta}{\sin\varphi\,\sin^2\psi}\,(\mathcal{R}+i\mathcal{I}).
$$

### 3.3 精确投影条件

相位方程（模 $\pi$ 符号由中心元自动吸收）：

$$
\frac{\mathcal{I}}{\mathcal{R}}=\tan\frac{3\pi}{n}.
$$

模长方程：

$$
\frac{2|S|^3|\cos\beta|}{|\sin\varphi|\sin^2\psi}\sqrt{\mathcal{R}^2+\mathcal{I}^2}=\frac{\cos(\pi/n)}{4}.
$$

---

## 4. 固定锚点的解析闭式解

### 4.1 锚点选取

取固定几何参数

$$
\boxed{\beta=\frac{\pi}{3}},\qquad\boxed{\psi=\frac{\pi}{4}}.
$$

此时投影公式简化为

$$
\mathcal{R}=\frac{3}{4}\cos\varphi-\frac{1}{2}\sin\varphi,\qquad
\mathcal{I}=-\frac{\sqrt{3}}{4}\cos\varphi.
$$

### 4.2 相位方程的显式解

令 $T=\tan(3\pi/n)$，由 $\mathcal{I}/\mathcal{R}=T$ 解得

$$
\boxed{\tan\varphi=\frac{3T+\sqrt{3}}{2T}}.
$$

对 $n=6$，$T\to\infty$，取极限得 $\tan\varphi=3/2$，与已知解析解一致。

### 4.3 模长方程的显式解

由模长条件解得

$$
\boxed{|S|=\sqrt[3]{\frac{|\sin\varphi|\cos(\pi/n)}{8\sqrt{\mathcal{R}^2+\mathcal{I}^2}}}}.
$$

$S$ 的符号由 $\pi(P)=$ 精确目标值确定（三种立方根中取使相位匹配的实根）。

### 4.4 生成元系数

$$
\boxed{
\begin{aligned}
a_2&=a_3=\frac{S}{\sin\psi}=S\sqrt{2},\\[4pt]
a_1&=\frac{-2S\cos\beta}{\sin\varphi}=-\frac{S}{\sin\varphi},\\[4pt]
a_0&=1-a_1\cos\varphi-2a_2\cos\psi.
\end{aligned}
}
$$

---

## 5. $n=6$ 的闭合解析构造

### 5.1 参数

取 $n=6$，$g=\gcd(3,6)=3$，$n/g=2$。二维目标为

$$
z=-\frac{\cos(\pi/6)}{4}\,e^{i\pi/2}=-\frac{\sqrt{3}}{8}\,i,\qquad
Z_1=z^2=-\frac{3}{64}.
$$

选取

$$
\beta=\frac{\pi}{3},\quad\psi=\frac{\pi}{4},\quad\varphi=\arctan\frac{3}{2},
$$

即 $\cos\varphi=2/\sqrt{13}$，$\sin\varphi=3/\sqrt{13}$。

### 5.2 投影计算

代入得 $\mathcal{R}=0$，$\mathcal{I}=-\sqrt{3}/(2\sqrt{13})$，故

$$
\pi(Q_1Q_2Q_3)=-i\frac{\sqrt{3}}{2\sqrt{13}}.
$$

环流投影

$$
\pi(P)=\frac{-S^3}{(3/\sqrt{13})(1/2)}\left(-i\frac{\sqrt{3}}{2\sqrt{13}}\right)=i\frac{\sqrt{3}}{3}S^3.
$$

令 $\pi(P)=-i\sqrt{3}/8$，解得

$$
\boxed{S=-\frac{\sqrt[3]{3}}{2}}.
$$

### 5.3 生成元系数

$$
\boxed{
\begin{aligned}
a_2&=a_3=-\frac{\sqrt[3]{3}}{\sqrt{2}},\\[4pt]
a_1&=\frac{\sqrt[3]{3}\sqrt{13}}{6},\\[4pt]
a_0&=1+\frac{2\sqrt[3]{3}}{3}.
\end{aligned}
}
$$

### 5.4 公理验证

**(A1)** 矢部求和：

$$
a_1\sin\varphi\,\mathbf{u}_1+a_2\sin\psi\,\mathbf{u}_2+a_3\sin\psi\,\mathbf{u}_3
=\frac{\sqrt[3]{3}}{2}j-\frac{\sqrt[3]{3}}{4}j-\frac{\sqrt[3]{3}}{4}k
+\left(-\frac{\sqrt[3]{3}}{4}j+\frac{\sqrt[3]{3}}{4}k\right)=\mathbf{0}.
$$

实部求和：

$$
a_0+a_1\cos\varphi+2a_2\cos\psi
=\left(1+\frac{2\sqrt[3]{3}}{3}\right)+\frac{\sqrt[3]{3}}{3}-\sqrt[3]{3}=1.
$$

**(A2)** $\pi(P)=-i\sqrt{3}/8$，与二维目标完全一致。

**(A3)** $\pi(P^2)=(-i\sqrt{3}/8)^2=-3/64=Z_1$，中心元精确匹配。

---

## 6. 主定理：DCG 维度提升

**定理（DCG 维度提升）**  
设二维 DCG 系统由整数 $n\ge 3$ 参数化，其算子族为
$$
M_1=\frac12 I,\qquad M_2=-\frac12 R_n,\qquad M_3=\frac12(I+R_n),
$$
其中 $R_n\in\mathrm{SO}(2)$ 为旋转 $2\pi/n$ 的矩阵。记环流算子
$$
P=M_1M_2M_3=-\frac{\cos(\pi/n)}{4}\,e^{i3\pi/n},
$$
并令 $g=\gcd(3,n)$，则本原中心元
$$
Z_1=P^{\,n/g}=(-1)^{\frac{n+3}{g}}\Bigl(\frac{\cos(\pi/n)}{4}\Bigr)^{\!n/g}\,I.
$$

则存在一族三维四元数算子 $\{\tilde M_0,\tilde M_1,\tilde M_2,\tilde M_3\}\subset\mathbb{H}$，满足：

1. **单位分解**：$\displaystyle\sum_{k=0}^3\tilde M_k=I$；
2. **环流投影匹配**：$\pi(\tilde M_1\tilde M_2\tilde M_3)=P$；
3. **中心元保持**：$\pi(\tilde Z_1)=Z_1$，其中 $\tilde Z_1=(\tilde M_1\tilde M_2\tilde M_3)^{n/g}$，

而投影映射 $\pi:\mathbb{H}\to\mathbb{C}$ 定义为 $\pi(a+bi+cj+dk)=a+bi$。

**唯一性（刚性）**  
在标准轴对称配置
$$
\mathbf{u}_1=j,\quad
\mathbf{u}_2=\cos\beta\,j+\sin\beta\,k,\quad
\mathbf{u}_3=\cos\beta\,j-\sin\beta\,k,
$$
并固定几何锚点
$$
\beta=\frac{\pi}{3},\qquad \psi=\frac{\pi}{4},
$$
则上述提升对每个 $n\ge 3$ 唯一确定。此时 $\tilde M_k=a_k Q_k$（$Q_k=\cos\phi_k+\mathbf{u}_k\sin\phi_k$），其中

- 令 $T=\tan\dfrac{3\pi}{n}$，则
  $$
  \tan\varphi=\frac{3T+\sqrt{3}}{2T},\qquad \phi_2=\phi_3=\psi=\frac{\pi}{4},\qquad \phi_1=\varphi;
  $$

- 记 $S=a_2\sin\psi=a_3\sin\psi$，其模长为
  $$
  |S|=\sqrt[3]{\frac{|\sin\varphi|\cos(\pi/n)}{8\sqrt{\mathcal{R}^2+\mathcal{I}^2}}},
  $$
  其中
  $$
  \mathcal{R}=\frac34\cos\varphi-\frac12\sin\varphi,\qquad
  \mathcal{I}=-\frac{\sqrt3}{4}\cos\varphi,
  $$
  而 $S$ 的符号由 $\pi(\tilde M_1\tilde M_2\tilde M_3)=P$ 的相位唯一锁定；

- 生成元系数由初等代数数显式给出：
  $$
  a_2=a_3=\frac{S}{\sin\psi}=S\sqrt2,\qquad
  a_1=-\frac{S}{\sin\varphi},\qquad
  a_0=1-a_1\cos\varphi-2a_2\cos\psi.
  $$

**注**  
条件 (3) 并非独立公理：由 $\pi(\tilde P)=P$ 及 $\pi$ 为同态，归纳可得 $\pi(\tilde P^k)=P^k$ 对所有正整数 $k$ 成立，故中心元保持是环流投影匹配的必然推论。

---

**定义（三进制能谱母哈密顿量）**  
在 Hilbert 空间 $\mathcal{H}=\ell^2(\mathbb{N}_{\ge 3})$ 上，取标准正交基 $\{|n\rangle\}_{n=3}^\infty$。定义三进制相位算子
$$
\Theta|n\rangle=\omega^n|n\rangle,\qquad \omega=e^{2\pi i/3},
$$
满足 $\Theta^3=I$，其本征值 $1,\omega,\omega^2$ 分别对应 $n\equiv 0,1,2\pmod 3$。

定义三格点平移算子
$$
\hat{T}_3|n\rangle=|n+3\rangle,
$$
显然 $[\hat{T}_3,\Theta]=0$。

**DCG 母哈密顿量** $H_{\mathrm{DCG}}:\mathcal{H}\to\mathcal{H}$ 定义为
$$
\boxed{H_{\mathrm{DCG}}=H_0+\gamma H_1},\qquad \gamma>0,
$$
其中
$$
\begin{aligned}
H_0&=\sum_{n=3}^{\infty}E_n\,|n\rangle\langle n|,\\[4pt]
E_n&=\frac{n}{g}\ln\!\left(\frac{4}{\cos(\pi/n)}\right),\qquad g=\gcd(3,n),\\[6pt]
H_1&=\hat{T}_3+\hat{T}_3^{\dagger}
    =\sum_{n=3}^{\infty}\bigl(|n\rangle\langle n+3|+|n+3\rangle\langle n|\bigr).
\end{aligned}
$$

**定理（三进制超选择定则与能带结构）**  
$H_{\mathrm{DCG}}$ 与 $\Theta$ 对易：
$$
[H_{\mathrm{DCG}},\Theta]=0.
$$

因此 $\mathcal{H}$ 分解为三个相互正交的 sector：
$$
\mathcal{H}=\mathcal{H}_0\oplus\mathcal{H}_1\oplus\mathcal{H}_2,
\qquad \mathcal{H}_r=\operatorname{span}\{|n\rangle:n\equiv r\pmod 3\}.
$$

在每个 sector 内，$H_{\mathrm{DCG}}$ 的本征能谱形成独立的三进制能带。大 $n$ 时 $E_n\approx \frac{n}{g}\ln 4$，故：

- **Sector $r=0$**（$n=3m$，$g=3$）：能级基线最低，平均间距 $\Delta E_0\approx\ln 4$，态密度为其余两 sector 的三倍。
- **Sector $r=1,2$**（$n=3m+1,3m+2$，$g=1$）：能级基线较高，平均间距 $\Delta E_{1,2}\approx 3\ln 4$，两 sector 仅由相位 $\omega$ 与 $\omega^2$ 区分。

**注**  
能量本征值 $E_n$ 由二维 DCG 本原中心元模长 $-\ln|Z_1^{(n)}|$ 直接读出，故该哈密顿量的谱结构完全由离散共形几何的代数数据锁定，不含任何自由拟合参数。