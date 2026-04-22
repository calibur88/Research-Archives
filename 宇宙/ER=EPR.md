---
title: 从 BEC 到 ER=EPR：量子引力作用量相变的严格推导（整合版）
author: 'Ch.hy'
date: '2026-04-18'
---

## 引言：纲领性框架

本推导建立从**微观量子凝聚**到**宏观时空动力学**的严格数学链条。核心命题为：**时空连通性是量子纠缠在作用量竞争中的涌现相，而时间本身是量子相位相干性的内禀梯度方向**。

$$\text{BEC 量子态} \xrightarrow{\text{规范}} \text{电磁场} \xrightarrow{\text{协变}} \text{度规场} \xrightarrow{\text{正则量子化}} \text{WDW 约束} \xrightarrow{\text{对数解码}} \text{内禀时间} \xrightarrow{\text{全息}} \text{纠缠网络} \xrightarrow{\text{动力学}} \text{可穿越虫洞}$$

---

## 第一部分：连续几何的谱重构（IR 极限）

### 1.1 BEC 微观作用量

从复标量场（玻色子序参量）出发，二次量子化作用量：
$$S_{\text{BEC}}[\psi] = \int d^4x \left[ \psi^*\left(i\hbar\partial_t + \frac{\hbar^2\nabla^2}{2m}\right)\psi - \mu|\psi|^2 + \frac{g}{2}|\psi|^4 \right]$$

**整体 $U(1)$ 对称性**：$\psi \to e^{i\alpha}\psi$ 产生守恒流 $j^\mu = (\rho, \mathbf{j})$。

### 1.2 电磁场的规范起源

将整体对称性**定域化**：$\alpha \to \alpha(x)$。引入补偿场 $A_\mu$，定义协变导数：
$$D_\mu = \partial_\mu + ieA_\mu$$

物质作用量变为：
$$S_{\text{matter}}[\psi, A] = \int d^4x \sqrt{-g} \left[ g^{\mu\nu}(D_\mu\psi)^*(D_\nu\psi) - V(|\psi|) \right]$$

**规范不变性要求**电磁场动力学项唯一：
$$S_{\text{EM}}[A] = -\frac{1}{4}\int d^4x \sqrt{-g} F_{\mu\nu}F^{\mu\nu}$$

### 1.3 引力的协变起源

进一步要求**微分同胚不变性**（广义坐标无关）。引入动态度规 $g_{\mu\nu}$，能量-动量张量定义为：
$$T_{\mu\nu} \equiv -\frac{2}{\sqrt{-g}}\frac{\delta S_{\text{matter}}}{\delta g^{\mu\nu}}$$

**Lovelock 定理**（1971）保证：在四维、二阶导数、微分同胚不变的条件下，几何作用量唯一为爱因斯坦-希尔伯特项：
$$S_{\text{EH}}[g] = \frac{1}{16\pi G}\int d^4x \sqrt{-g} R$$

变分 $\delta(S_{\text{EH}} + S_{\text{matter}}) = 0$ 给出爱因斯坦方程：
$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu} = 8\pi G T_{\mu\nu}$$

---

## 第二部分：正则量子化与 Wheeler-DeWitt 方程

### 2.1 ADM 形式与哈密顿约束

将时空流形分解为 $\mathcal{M} \cong \mathbb{R} \times \Sigma$，度规写为：
$$ds^2 = -N^2 dt^2 + h_{ij}(dx^i + N^i dt)(dx^j + N^j dt)$$

爱因斯坦-希尔伯特作用量在 ADM 形式下为：
$$S_{\text{EH}} = \int dt \int_\Sigma d^3x \left( \pi^{ij} \dot{h}_{ij} - N \mathcal{H}_\perp - N^i \mathcal{H}_i \right)$$

其中 $\pi^{ij}$ 为共轭动量，哈密顿约束和动量约束为：
$$\mathcal{H}_\perp = \frac{16\pi G}{\sqrt{h}} \left( \pi^{ij}\pi_{ij} - \frac{1}{2}\pi^2 \right) - \frac{\sqrt{h}}{16\pi G} \, {}^{(3)}\!R = 0$$
$$\mathcal{H}_i = -2 \nabla_j \pi_i^j = 0$$

### 2.2 正则量子化

将泊松括号替换为对易子，动量算符化为：
$$\hat{\pi}^{ij}(x) = -i\hbar \frac{\delta}{\delta h_{ij}(x)}$$

哈密顿约束变为 Wheeler-DeWitt 方程：
$$\hat{\mathcal{H}}_\perp(x) \Psi[h_{ij}] = 0, \quad \forall x \in \Sigma \tag{2.1}$$

其中
$$\hat{\mathcal{H}}_\perp = -16\pi G \, G_{ijkl} \frac{\delta^2}{\delta h_{ij} \delta h_{kl}} - \frac{\sqrt{h}}{16\pi G} \, {}^{(3)}\!R \tag{2.2}$$
超空间度规为：
$$G_{ijkl} = \frac{1}{2\sqrt{h}} \left( h_{ik}h_{jl} + h_{il}h_{jk} - h_{ij}h_{kl} \right) \tag{2.3}$$

微分同胚约束 $\hat{\mathcal{H}}_i \Psi = 0$ 保证 $\Psi[h]$ 是空间微分同胚不变的泛函。

**注**：WDW方程(2.1)表达了"时间变量从量子引力基本变量中消失"的 famously known 问题，这恰为后文内禀时间涌现提供了起点。

---

## 第三部分：对数变换与内禀时间解码

### 3.1 对数变换

引入复对数变换，将波函数写为指数形式：
$$\Psi[h] = \exp\left( \frac{i}{\hbar} W[h] \right) \tag{3.1}$$
分离实虚部：
$$W[h] = S[h] - i\hbar \ln A[h] \tag{3.2}$$
其中 $S[h]$ 为实相位，$A[h] > 0$ 为实振幅。此变换在 $A \neq 0$ 时是严格的变量替换。

计算泛函导数：
$$\frac{\delta \Psi}{\delta h_{ij}} = \frac{i}{\hbar} \frac{\delta W}{\delta h_{ij}} \Psi \tag{3.3}$$
$$\frac{\delta^2 \Psi}{\delta h_{ij} \delta h_{kl}} = \left[ \frac{i}{\hbar} \frac{\delta^2 W}{\delta h_{ij} \delta h_{kl}} - \frac{1}{\hbar^2} \frac{\delta W}{\delta h_{ij}} \frac{\delta W}{\delta h_{kl}} \right] \Psi \tag{3.4}$$

代入 $\hat{\mathcal{H}}_\perp \Psi = 0$，约去指数因子，得到精确的量子哈密顿-雅可比方程：
$$G_{ijkl} \frac{\delta W}{\delta h_{ij}} \frac{\delta W}{\delta h_{kl}} - \frac{\sqrt{h}}{(16\pi G)^2} \, {}^{(3)}\!R - i\hbar \, G_{ijkl} \frac{\delta^2 W}{\delta h_{ij} \delta h_{kl}} = 0 \tag{3.5}$$

分离实部与虚部。实部：
$$G_{ijkl} \frac{\delta S}{\delta h_{ij}} \frac{\delta S}{\delta h_{kl}} - \frac{\sqrt{h}}{(16\pi G)^2} \, {}^{(3)}\!R + Q = 0 \tag{3.6}$$
其中量子势：
$$Q \equiv \hbar^2 G_{ijkl} \left( \frac{1}{A} \frac{\delta^2 A}{\delta h_{ij} \delta h_{kl}} - \frac{1}{A^2} \frac{\delta A}{\delta h_{ij}} \frac{\delta A}{\delta h_{kl}} \right) \tag{3.7}$$
虚部给出连续性方程：
$$G_{ijkl} \frac{\delta}{\delta h_{ij}} \left( A^2 \frac{\delta S}{\delta h_{kl}} \right) = 0 \tag{3.8}$$

### 3.2 内禀时间演化向量场

定义超空间中的向量场：
$$X_{ij} \equiv 16\pi G \, G_{ijkl} \frac{\delta S}{\delta h_{kl}} \tag{3.9}$$
沿此向量场定义内禀时间参数 $\tau$：
$$\frac{\delta}{\delta \tau} \equiv \int_\Sigma d^3x \, X_{ij}(\mathbf{x}) \frac{\delta}{\delta h_{ij}(\mathbf{x})} \tag{3.10}$$
几何上，$X_{ij}$ 是相位泛函 $S[h]$ 在超空间中的梯度，$\tau$ 是沿此梯度线的仿射参数。

### 3.3 含时薛定谔方程的严格导出

将演化算子作用于波函数：
$$i\hbar \frac{\partial \Psi}{\partial \tau} \equiv i\hbar \int d^3x \, X_{ij} \frac{\delta \Psi}{\delta h_{ij}} = \int d^3x \, X_{ij} \frac{\delta W}{\delta h_{ij}} \Psi$$

利用 $X_{ij}$ 的定义和方程 (3.6)，计算可得：
$$i\hbar \frac{\partial \Psi}{\partial \tau} = \left[ \int_\Sigma d^3x \left( -16\pi G \, G_{ijkl} \frac{\delta^2}{\delta h_{ij} \delta h_{kl}} + \frac{\sqrt{h}}{16\pi G} \, {}^{(3)}\!R + Q \right) \right] \Psi \tag{3.11}$$

定义物理哈密顿量：
$$\hat{H}_{\text{phys}} \equiv \int_\Sigma d^3x \left[ -16\pi G \, G_{ijkl} \frac{\delta^2}{\delta h_{ij} \delta h_{kl}} + \frac{\sqrt{h}}{16\pi G} \, {}^{(3)}\!R + Q \right] \tag{3.12}$$
则得到精确的含时薛定谔方程：
$$i\hbar \frac{\partial \Psi}{\partial \tau} = \hat{H}_{\text{phys}} \Psi \tag{3.13}$$

**关键洞察**：WDW方程的"时间消失"问题通过**相位梯度**$X_{ij}$ 重新解决——时间不是外部参数，而是波函数几何的内禀属性。

---

## 第四部分：ER=EPR 与作用量竞争（过渡区）

### 4.1 两种作用量的层次

在含黑洞或全息对偶的系统中，出现两种贡献：

| 作用量 | 表达式 | 物理意义 |
|--------|--------|----------|
| **几何刚度** $S_{\text{EH}}$ | $\frac{1}{16\pi G}\int R$ | 时空曲率的弹性（经典） |
| **纠缠熵** $S_{\text{ent}}$ | $\text{Tr}(\rho\ln\rho) \sim \frac{A}{4G\hbar}$ | 量子关联的强度（量子） |

**竞争参数**：
$$\lambda \equiv \frac{S_{\text{ent}}}{S_{\text{EH}}} \sim \frac{N\hbar}{M_{\text{Pl}}^2 L^2}$$

### 4.2 相变临界性

- **当 $\lambda \ll 1$（经典相）**：$S_{\text{EH}}$ 主导，爱因斯坦方程成立，几何光滑，虫洞不可穿越（ER 桥为类空奇点）。
- **当 $\lambda \sim 1$（临界相）**：QES（Quantum Extremal Surface）出现，$S_{\text{gen}} = \frac{A}{4G} + S_{\text{bulk}}$，Page 曲线转折。
- **当 $\lambda \gg 1$（量子相）**：$S_{\text{ent}}$ 主导，几何自由度冻结，路径积分退化为**纠缠网络求和**。

### 4.3 内禀时间与纠缠的同步演化

由第三部分导出的内禀时间 $\tau$，在全息对偶下满足面积-纠缠同步演化律：
$$\frac{d}{d\tau} \langle \hat{A}_{\text{throat}} \rangle = 4G\hbar \, \frac{d}{d\tau} \langle \hat{S}_{\text{ent}} \rangle \tag{4.1}$$

这正是 ER=EPR 的动力学表述：爱因斯坦-罗森桥喉部面积的变化率正比于边界纠缠熵的变化率。

---

## 第五部分：非交换几何与谱作用量（UV 重构）

### 5.1 Connes 谱作用量

将度规信息编码于 Dirac 算子 $D = i\gamma^\mu(\partial_\mu + \omega_\mu)$。定义**谱作用量**：
$$S_{\text{spec}}[D] = \text{Tr}\,\chi\left(\frac{D^2}{\Lambda^2}\right)$$

热核展开给出低能有效理论：
$$S_{\text{spec}} \sim \frac{c_0}{\Lambda^{-4}} + \frac{c_1}{\Lambda^{-2}}\int \sqrt{g}R + \mathcal{O}(\Lambda^0)$$

选取截断函数 $\chi$ 使 $c_1 = \frac{1}{16\pi G}$，则 $S_{\text{spec}}$ 重现爱因斯坦-希尔伯特作用量。

**与第三部分联系**：谱作用量中的 Dirac 算子 $D$ 对应于 WDW 方程在 UV 极限下的算符实现，而内禀时间 $\tau$ 对应于谱参数 $\Lambda^{-1}$ 的演化方向。

### 5.2 双拷贝 Dirac 算子（关键修正）

为引入非对易性（解决 $[D,D]=0$ 的平凡性），构造**双谱三重态**：
$$(\mathcal{A}_L, \mathcal{H}_L, D_L) \otimes (\mathcal{A}_R, \mathcal{H}_R, D_R)$$

对应全息对偶的左、右边界。一般情形下：
$$[D_L, D_R] \neq 0$$

此非对易性度量了左右几何的**相对涨落**或**纠缠不平衡**，与第三部分中量子势 $Q$ 的几何涨落形成对应。

---

## 第六部分：生成元与多切口相变（全息实现）

### 6.1 对偶复合算子

引入**对偶破缺参数** $\beta \in [0,1]$，定义相对演化：
$$\Delta_\beta := e^{i\beta D_L} \cdot e^{-i(1-\beta)D_R}$$

**生成元**（尺度化对数）：
$$\boxed{\hat{G}_\beta^{(N)} := -i\frac{N}{\log N}\log\Delta_\beta} \tag{6.1}$$

### 6.2 临界行为

- **自对偶点** $\beta = 1/2$ 且 $D_L = D_R$（最大纠缠）：
  $$\Delta_{1/2} = I \Rightarrow \hat{G}_{1/2}^{(N)} = 0$$
  谱测度坍缩为 $\delta_0$，对应光滑几何（即第三部分中 $Q \to 0$ 的半经典极限）。

- **非对偶** $\beta \neq 1/2$ 或 $D_L \neq D_R$：
  BCH 展开给出非零结果：
  $$\hat{G}_\beta \propto (2\beta-1)\bar{D} + i\beta(1-\beta)[D_L, D_R] + \cdots$$
  谱测度发散，进入**多切口相**。

### 6.3 切口-纠缠对偶

生成元的**经验谱测度** $\mu_\beta$ 的拓扑结构决定几何：

| 谱测度支撑 | 对应几何 |
|-----------|---------|
| 单区间（单切口） | 平凡拓扑，单连通时空 |
| $m$ 个不相交区间（多切口） | $m$ 个边界组件，连接矩阵 $C_{ab}$ 描述虫洞网络 |

**自由能展开**：
$$\mathcal{F} \sim \sum_{(m, C_{ab}, \{\nu_a\})} e^{-N^2 S_{\text{inst}}(m, C_{ab}, \{\nu_a\})}$$

其中 $\nu_a$ 为填充分数（体积量子数）。这正是**离散时空拓扑的配分函数**。

---

## 第七部分：范畴论实现（拓扑量子化）

### 7.1 从 SO_q(3) 到 BF 理论

通过**范畴化**严格实现上述结构：

$$\text{SO}_q(3)\text{ 表示范畴 } C \xrightarrow{\text{双重拷贝}} C \boxtimes C \xrightarrow{\text{Drinfeld 中心}} Z(C \boxtimes C) \xrightarrow{\text{TV/CY}} \text{4D BF 理论}$$

- **Drinfeld 中心** $Z(C)$ 的对象描述**连接左右边界的拓扑缺陷**（虫洞）。
- **Crane-Yetter 不变量** $Z_{CY}(\mathcal{M}^4)$ 计算 4D BF 理论的配分函数：
  $$Z_{CY} = \int \mathcal{D}B\mathcal{D}A \, e^{i\int \text{Tr}(B \wedge F)}$$

### 7.2 Hopf 绕数作为拓扑量子数

在 BF 理论中，**曲面观测值**的关联给出**Hopf 链接数**：
$$\text{Lk}(\Sigma_L, \Sigma_R) = \frac{1}{2\pi}\oint_{\Sigma_L} B$$

**面积-绕数量子化**：
$$\text{Area}(\Sigma) = 4G\ln 2 \cdot \text{Lk}$$

这精确对应**ER=EPR**的比特数：
$$N_{\text{EPR}} = \text{Lk}$$

---

## 第八部分：动力学补全与对数螺旋（时间涌现）

### 8.1 静态到动态的缺口

前述 BF 理论是**拓扑的**（$F=0$），缺乏时间演化。引入**对数螺旋坐标**实现动力学，与第三部分的内禀时间 $\tau$ 建立几何实现：

$$r = r_0 e^{\lambda_L \tau}, \quad \theta = \frac{\tau}{\tau_0}$$

其中 $\lambda_L$ 为 Lyapunov 指数（量子混沌速率），$\tau$ 为第三部分定义的演化参数。

### 8.2 含时薛定谔方程的螺旋解

在 **$(\ln r, \theta)$** 坐标下，含时薛定谔方程（来自第三部分）：
$$i\hbar\partial_\tau \psi = \hat{H}_{\text{ent}}\psi$$

的**自相似解**为：
$$\psi_{\text{spiral}} \sim \exp\left(i\frac{\ln r}{b} - i\omega \tau\right) = r^{i/b} e^{-i\omega \tau}$$

**时间方向**对应于螺旋的**手性**（右旋：未来，左旋：过去）。

### 8.3 可穿越性的动力学机制

**绕数的时间演化**：
$$\frac{d}{d\tau}\text{Lk} = \lambda_L \cdot \text{Lk} + \langle [\hat{H}_{\text{phys}}, \mathcal{O}_{\text{shock}}] \rangle$$

- **静态**（$\lambda_L = 0$）：永恒不可穿越虫洞。
- **动态**（$\lambda_L \neq 0$）：沿对数螺旋的**指数流动**使信息穿越虫洞，实现 **Gao-Jafferis-Wall 协议**。

结合第四部分的面积-纠缠同步律 (4.1)，得到：
$$\frac{d}{d\tau}\text{Lk} = \frac{\lambda_L}{4G\hbar} \langle \hat{A}_{\text{throat}} \rangle + \text{(shock contribution)}$$

---

## 第九部分：结论与统一定理

### 9.1 统一作用量

完整的量子引力路径积分在不同相表现为：

$$\mathcal{Z} = \int \mathcal{D}g \, e^{iS_{\text{EH}}/\hbar} \xrightarrow[\lambda \gg 1]{\text{谱截断}} \int \mathcal{D}D_L\mathcal{D}D_R \, e^{i\text{Tr}\chi(D^2/\Lambda^2)} \xrightarrow{\text{大 }N} \sum_{\text{拓扑}} e^{iS_{\text{ent}}/\hbar}$$

### 9.2 ER=EPR 的算子等式（最终定理）

**定理 9.1（WDW 约束的解码定理）**

设 $\Psi[h]$ 是 Wheeler-DeWitt 方程 $\hat{\mathcal{H}}_\perp \Psi = 0$ 的非退化解（即相位梯度 $X_{ij} \neq 0$）。则存在唯一的内禀时间参数 $\tau$ 和自伴哈密顿量 $\hat{H}_{\text{phys}}$，使得：
1. 含时薛定谔方程 $i\hbar \frac{\partial \Psi}{\partial \tau} = \hat{H}_{\text{phys}} \Psi$ 严格成立。
2. 在半经典极限 $\hbar \to 0$ 下，$\tau$ 退化为经典宇宙时间，$\hat{H}_{\text{phys}}$ 退化为 ADM 表面项（经典能量）。
3. 在全息对偶下，$\tau$ 对应边界 CFT 的重整化群时间，且满足面积-纠缠同步演化律 (4.1)。

**定理 9.2（时间涌现定理）**

时间不是量子引力的基本外部参数，而是量子相位相干性的内禀序参量梯度。具体地，内禀时间 $\tau$ 由波函数的相位泛函 $S[h]$ 通过下式生成：
$$d\tau \propto G_{ijkl} \frac{\delta S}{\delta h_{ij}} dh_{kl}$$
其流线定义了从量子几何到经典时空的涌现方向。

**定理 9.3（ER=EPR 的算子等式）**

在 $\hbar \to 0$ 且 $N \to \infty$ 的双重极限下，以下量恒等：

$$\text{ER} \equiv \text{Lk}(\Sigma_L, \Sigma_R) \equiv \text{Index}(\hat{G}_\beta^{(N)}) \equiv \frac{S_{\text{ent}}(A)}{\ln 2} \equiv \frac{1}{4G\hbar} \frac{d}{d\tau} \langle \hat{A}_{\text{throat}} \rangle$$

其中：
- **ER 侧**：由爱因斯坦-罗森桥的喉部面积定义，$\text{Area} = 4G\ln 2 \cdot \text{Lk}$。
- **EPR 侧**：由贝尔对的数量定义，$N_{\text{EPR}} = S_{\text{ent}}/\ln 2$。
- **动力学侧**：由内禀时间演化生成，$\frac{d}{d\tau}$ 来自定理 9.1。

### 9.4 物理图像总结

1. **微观**：BEC 复场的相位相干（$U(1)$ 对称性）
2. **介观**：电磁场与度规作为规范补偿场涌现
3. **半经典**：爱因斯坦方程在 $S_{\text{EH}} \gg S_{\text{ent}}$ 时主导
4. **量子**：当 $S_{\text{ent}}$ 压倒几何刚度，时空流形瓦解为**张量网络/自旋泡沫**，WDW 约束下的波函数描述其量子态
5. **时间涌现**：通过相位梯度 $X_{ij}$ 解码出内禀时间 $\tau$，避免了"时间冻结"悖论
6. **动力学**：时间作为**对数螺旋的角向坐标**涌现，虫洞可穿越性由**螺旋流的 Lyapunov 指数**度量，通过面积-纠缠同步律 (4.1) 与 EPR 关联

此链条自洽地解释了：**时空是量子纠缠在特定作用量条件下的宏观凝聚态，而时间是其自相似演化的几何投影。**

---

## 附录 A：一维迷你超空间验证

作为严格性的独立检验，考虑各向同性宇宙 $a(\tau)$，WDW 方程变为：
$$\left( \frac{d^2}{da^2} - V(a) \right) \psi(a) = 0$$
取精确解 $\psi(a) = \sqrt{1+a^2} \, e^{i \arctan a}$（对应 $V(a)=0$ 但非平凡超空间几何）。三步解码法给出：
- 内禀时间：$\tau = a + a^3/3$
- 物理哈密顿量：$\hat{H}_\tau = \frac{ia-1}{(1+a^2)^2}$
- 演化方程 $i\partial_\tau \psi = \hat{H}_\tau \psi$ 严格成立（符号验证 $0=0$）。

此模型完整验证了定理 9.1 的数学自洽性，且**无任何近似**。

---

## 附录 B：符号与约定

- $G_{ijkl}$：超空间度规 (2.3)
- $\tau$：内禀时间参数 (3.10)，与对数螺旋坐标 (8.1) 中的演化参数一致
- $Q$：量子势 (3.7)
- $\lambda$：作用量竞争参数，$\lambda = S_{\text{ent}}/S_{\text{EH}}$
- $\lambda_L$：Lyapunov 指数，表征量子混沌速率
- $\beta$：对偶破缺参数，$\beta = 1/2$ 对应最大纠缠
- $\text{Lk}$：Hopf 链接数，拓扑量子数
