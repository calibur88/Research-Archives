# 三维 DCG 的精确投影条件

## 摘要

设 $n\ge 3$ 为整数，$\zeta=e^{2\pi i/n}$。二维 DCG 代数由旋转矩阵生成，其环流算子 $P$ 的本原中心元 $Z_1=P^{n/\gcd(3,n)}$ 控制分形谱维数。本文以四元数 $\mathbb{H}$ 为载体，通过单位分解与中心元投影等式两条公理，构造了二维代数到三维非交换算子的提升。在锚点 $\beta=\pi/3,\psi=\pi/4$ 下给出唯一解析闭式解，并推广至任意 $\mathbb{Z}_m$ 对称性。投影映射 $\pi:\mathbb{H}\to\mathbb{C}$ 非乘法同态，但满足 $(\pi(\tilde P))^{n/g}=Z_1$，确保三维非交换信息无损编码于二维分形测度。进一步，在 Hilbert 空间 $\ell^2(\mathbb{N}_{\ge 3})$ 上建立 DCG 母哈密顿量 $H_{\mathrm{DCG}}$，其三进制超选择定则锁定能谱为三个独立 sector，基线能量由 $-\ln|Z_1^{(n)}|$ 直接读出，不含任何自由拟合参数。

---

## 1. 二维 DCG 代数

### 1.1 旋转矩阵与生成元

设整数 $n\ge 3$，记 $\zeta=e^{2\pi i/n}$。二维实旋转矩阵 $R\in\mathrm{SO}(2)$ 满足特征方程

$$R^2-2\cos\frac{2\pi}{n}R+I=0,$$

其特征值为 $\zeta,\bar\zeta$。定义三个线性算子

$$M_1=\frac12 I,\qquad M_2=-\frac12 R,\qquad M_3=\frac12(I+R),$$

它们两两对易，满足单位分解 $M_1+M_2+M_3=I$，生成的实代数

$$\mathcal{A}_n=\mathbb{R}[R]/(R^2-2\cos\tfrac{2\pi}{n}R+I)\cong\mathbb{C}.$$

### 1.2 算子字与标量退化

对非负整数 $a,b,c$，定义算子字 $W(a,b,c)=M_1^a M_2^b M_3^c$。在特征基下对角化，其特征值为

$$\lambda_\zeta=2^{-(a+b+c)}(-1)^b\zeta^b(1+\zeta)^c,\quad
\lambda_{\bar\zeta}=\overline{\lambda_\zeta}.$$

**定理 1.1（标量退化）** 算子字退化为标量矩阵 $W(a,b,c)\in\mathbb{R}I$ 当且仅当

$$2b+c\equiv 0\pmod n.$$

此时令 $2b+c=nk$，标量值为

$$W(a,b,c)=(-1)^{b+k}2^{-(a+b+c)}\bigl(2\cos\tfrac{\pi}{n}\bigr)^{nk-2b}\,I.$$

### 1.3 环流算子与中心元

定义环流算子 $P=M_1M_2M_3$。直接计算得

$$P=-\frac{\zeta(1+\zeta)}{8}=-\frac{\cos(\pi/n)}{4}\,e^{i3\pi/n}.$$

令 $g=\gcd(3,n)$，则 $P^{\,n/g}$ 的虚部归零，得到本原中心元

$$Z_1=P^{\,n/g}=(-1)^{\frac{n+3}{g}}\Bigl(\frac{\cos(\pi/n)}{4}\Bigr)^{\!n/g}\,I,$$

其模长 $|\lambda|=(\cos(\pi/n)/4)^{n/g}<1$。

### 1.4 分形维数

中心元的迭代生成自相似 IFS，每步分裂为 $3$ 个压缩比为 $|\lambda|$ 的副本，豪斯多夫维数为

$$d_f(n)=\frac{g\log 3}{n\log\bigl(4/\cos\frac{\pi}{n}\bigr)}.$$

当 $n\to\infty$ 时 $d_f(n)\to 0$。序列在极坐标 $(\theta_n,r_n)=(2\pi n/3,d_f(n))$ 下呈螺旋梯分布：$g=3$ 时 $\theta\equiv 0$，$g=1$ 时 $\theta$ 交替取 $2\pi/3,4\pi/3$。

---

## 2. 三维非交换提升

### 2.1 四元数生成元与投影

取四元数代数 $\mathbb{H}$，虚单位 $i,j,k$。设生成元

$$M_0=a_0 I,\qquad M_k=a_k Q_k\quad(k=1,2,3),$$

其中 $Q_k=\cos\phi_k+\mathbf{u}_k\sin\phi_k$ 为单位四元数，$\phi_k$ 为半角，轴 $\mathbf{u}_k\in\mathrm{span}\{j,k\}$ 为单位纯虚四元数。定义投影 $\pi:\mathbb{H}\to\mathbb{C}$ 为

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

其中 $\tilde P = \tilde M_1\tilde M_2\tilde M_3$，$g=\gcd(3,n)$，$Z_1$ 为二维本原中心元。此式要求先投影后取幂。当 $n/g=1$ 时退化为 $\pi(\tilde P)=P$。该等式完整确定所有 $\mathbb{Z}_2$ 规范自由度。

**注** (A2) 统合了环流投影匹配与中心元自洽两条要求。当 $n/g=1$（即 $n=3$）时退化为 $\pi(\tilde P)=P$；当 $n/g>1$ 时仅要求投影的幂次等于 $Z_1$，允许 $\pi(\tilde P)$ 与 $P$ 相差一个 $n/g$ 次单位根。该自由度由 $S$ 的三个实立方根分支覆盖，恰好锁定唯一解。

---

## 3. 对称假设与投影公式

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

---

## 4. 锚点下的解析解

### 4.1 锚点

取

$$\beta=\frac{\pi}{3},\qquad \psi=\frac{\pi}{4}.$$

此时

$$\mathcal{R}=\frac{3}{4}\cos\varphi-\frac{1}{2}\sin\varphi,\qquad
\mathcal{I}=-\frac{\sqrt{3}}{4}\cos\varphi.$$

### 4.2 参数确定

记 $T=\tan(3\pi/n)$。

**第一步** 由相位条件 $\mathcal{I}/\mathcal{R}=T$ 解得

$$\tan\varphi=\frac{3T+\sqrt{3}}{2T}.$$

**第二步** 由模长条件得

$$|S|=\sqrt[3]{\frac{|\sin\varphi|\cos(\pi/n)}{8\sqrt{\mathcal{R}^2+\mathcal{I}^2}}}.$$

**第三步** 利用中心元等式 $(\pi(\tilde P))^{n/g}=Z_1$ 从 $|S|$ 的三个立方根中唯一确定 $S$。

令 $S_k = |S| \cdot e^{i 2\pi k/3}$（$k=0,1,2$），对应 $\pi(\tilde P)_k = \pi(\tilde P)_0 \cdot e^{-i 2\pi k/3}$（因 $S^3$ 在分子）。则

$$(\pi(\tilde P)_k)^{n/g} = (\pi(\tilde P)_0)^{n/g} \cdot e^{-i 2\pi k (n/g)/3}.$$

由于 $\gcd(n/g,3)=1$，当 $k$ 遍历 $0,1,2$ 时指数遍历模 $3$ 完全剩余系，恰有一个 $k$ 使该因子抵消 $\pi(\tilde P)_0^{n/g}$ 与 $Z_1$ 之间的相位差，从而唯一确定 $S$ 的符号。

### 4.3 生成元系数

$$\begin{aligned}
a_2&=a_3=S\sqrt{2},\\
a_1&=-\frac{2S\cos\beta}{\sin\varphi}=-\frac{S}{\sin\varphi},\\
a_0&=1-a_1\cos\varphi-2a_2\cos\psi.
\end{aligned}$$

---

## 5. $n=6$ 的构造

$n=6$，$g=3$，$n/g=2$。二维目标

$$Z_1=\left(-\frac{\sqrt{3}}{8}i\right)^2=-\frac{3}{64}.$$

由 $\tan\varphi=3/2$ 得 $\cos\varphi=2/\sqrt{13},\sin\varphi=3/\sqrt{13}$。

计算得 $\mathcal{R}=0,\mathcal{I}=-\sqrt{3}/(2\sqrt{13})$，

$$\pi(\tilde P)=\frac{-S^3}{\sin\varphi}\mathcal{I}=i\frac{\sqrt{3}}{3}S^3.$$

代入 $(\pi(\tilde P))^2=Z_1$：

$$\left(i\frac{\sqrt{3}}{3}S^3\right)^2=-\frac{3}{64}\;\Longrightarrow\; S^6=\frac{9}{64}.$$

$S^3=\pm 3/8$。要求 $\pi(\tilde P)=P=-i\sqrt{3}/8$，故 $S^3=-3/8$，得

$$S=-\frac{\sqrt[3]{3}}{2}.$$

生成元系数：

$$\begin{aligned}
a_2&=a_3=-\frac{\sqrt[3]{3}}{\sqrt{2}},\\
a_1&=\frac{\sqrt[3]{3}\sqrt{13}}{6},\\
a_0&=1+\frac{2\sqrt[3]{3}}{3}.
\end{aligned}$$

验证公理：

**(A1)** 矢量部分：$a_1\sin\varphi\,\mathbf{u}_1+a_2\sin\psi(\mathbf{u}_2+\mathbf{u}_3)=\frac{\sqrt[3]{3}}{2}j-\frac{\sqrt[3]{3}}{2}j=\mathbf{0}$；

标量部分：$a_0+a_1\cos\varphi+2a_2\cos\psi = 1+\frac{2\sqrt[3]{3}}{3}+\frac{\sqrt[3]{3}}{3}-\sqrt[3]{3}=1$。

**(A2)** $\pi(\tilde P)=-i\sqrt{3}/8$，$(\pi(\tilde P))^2=-3/64=Z_1$。

---

## 6. 主定理

**定理 6.1（DCG 维度提升）** 对任意 $n\ge 3$，存在一族三维四元数算子 $\{\tilde M_0,\tilde M_1,\tilde M_2,\tilde M_3\}\subset\mathbb{H}$ 满足：

1. $\sum_{k=0}^3 \tilde M_k=I$；
2. $(\pi(\tilde M_1\tilde M_2\tilde M_3))^{n/g}=Z_1$，其中 $g=\gcd(3,n)$。

在轴对称配置 $\mathbf{u}_1=j,\mathbf{u}_2=\cos\beta\,j+\sin\beta\,k,\mathbf{u}_3=\cos\beta\,j-\sin\beta\,k$ 及锚点 $\beta=\pi/3,\psi=\pi/4$ 下，该提升唯一。投影 $\pi$ 非乘法同态，但满足 $(\pi(\tilde P))^{n/g}=Z_1$ 及 $\pi\circ\iota=\mathrm{id}$，其中 $\iota:\mathcal{A}_n\to\mathbb{H}$ 为二维环流算子到三维环流算子的自然嵌入 $P\mapsto\tilde P$。此嵌入是分裂满射 $\pi$ 的截面，保证三维非交换信息通过投影无损编码于二维分形测度。分形维数 $d_f(n)=\frac{g\log 3}{n\log(4/\cos(\pi/n))}$ 保持不变。

---

## 7. DCG 母哈密顿量与三进制超选择定则

### 7.1 母哈密顿量

在 Hilbert 空间 $\mathcal{H}=\ell^2(\mathbb{N}_{\ge 3})$ 上，取标准正交基 $\{|n\rangle\}_{n=3}^\infty$。定义三进制相位算子

$$\Theta|n\rangle=\omega^n|n\rangle,\qquad \omega=e^{2\pi i/3},$$

满足 $\Theta^3=I$，其本征值 $1,\omega,\omega^2$ 分别对应 $n\equiv 0,1,2\pmod 3$。

定义三格点平移算子

$$\hat{T}_3|n\rangle=|n+3\rangle,$$

满足 $[\hat{T}_3,\Theta]=0$。

DCG 母哈密顿量 $H_{\mathrm{DCG}}:\mathcal{H}\to\mathcal{H}$ 定义为

$$H_{\mathrm{DCG}}=H_0+\gamma H_1,\qquad \gamma>0,$$

其中

$$\begin{aligned}
H_0&=\sum_{n=3}^{\infty}E_n\,|n\rangle\langle n|,\\[4pt]
E_n&=\frac{n}{g}\ln\!\left(\frac{4}{\cos(\pi/n)}\right),\qquad g=\gcd(3,n),\\[6pt]
H_1&=\hat{T}_3+\hat{T}_3^{\dagger}
    =\sum_{n=3}^{\infty}\bigl(|n\rangle\langle n+3|+|n+3\rangle\langle n|\bigr).
\end{aligned}$$

### 7.2 三进制超选择定则

**定理 7.1（三进制超选择定则与能带结构）** $H_{\mathrm{DCG}}$ 与 $\Theta$ 对易：

$$[H_{\mathrm{DCG}},\Theta]=0.$$

因此 $\mathcal{H}$ 分解为三个相互正交的 sector：

$$\mathcal{H}=\mathcal{H}_0\oplus\mathcal{H}_1\oplus\mathcal{H}_2,
\qquad \mathcal{H}_r=\operatorname{span}\{|n\rangle:n\equiv r\pmod 3\}.$$

在每个 sector 内，$H_{\mathrm{DCG}}$ 的本征能谱形成独立的三进制能带。大 $n$ 时 $E_n\approx \frac{n}{g}\ln 4$，故：

- **Sector $r=0$**（$n=3m$，$g=3$）：能级基线最低，平均间距 $\Delta E_0\approx\ln 4$，态密度为其余两 sector 的三倍。
- **Sector $r=1,2$**（$n=3m+1,3m+2$，$g=1$）：能级基线较高，平均间距 $\Delta E_{1,2}\approx 3\ln 4$，两 sector 仅由相位 $\omega$ 与 $\omega^2$ 区分。

能量本征值 $E_n$ 由二维 DCG 本原中心元模长 $-\ln|Z_1^{(n)}|$ 直接读出，故该哈密顿量的谱结构完全由离散共形几何的代数数据锁定，不含任何自由拟合参数。

---

## 8. 广义 $\mathbb{Z}_m$ 对称性

令 $\beta=\pi/m\;(m\ge 3)$，$\psi=\pi/4$。投影公式为

$$\begin{aligned}
\mathcal{R}&=\cos\varphi\,\sin^2\beta-\sin\varphi\,\cos\beta,\\
\mathcal{I}&=-\cos\varphi\,\sin\beta\cos\beta.
\end{aligned}$$

相位方程：$\tan\varphi=\sin\beta(\tan\beta+1/T)$，$T=\tan(3\pi/n)$。

$|S|$ 由模长条件确定，$S$ 的精确值由 $(\pi(\tilde P))^{n/g}=Z_1$ 唯一锁定。生成元系数仍为

$$a_2=a_3=S\sqrt{2},\quad a_1=-\frac{2S\cos\beta}{\sin\varphi},\quad a_0=1-a_1\cos\varphi-\sqrt{2}a_2.$$

**定理 8.1（广义提升）** 对任意 $m\ge 3$ 与 $n\ge 3$，上述构造给出唯一三维提升，满足公理 (A1)(A2)。投影 $\pi$ 为分裂满射，分形维数与 $m$ 无关。