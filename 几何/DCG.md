# 离散构造几何（DCG）3.1

## 第一部分　声明与公理体系

### 1.1　声明

DCG 并非唯一几何学，而是特定合法性框架。公理 Ⅰ–Ⅴ 为选择性约束；凡未接受并满足者，DCG 不授予「几何闭包」之资格，亦不评价其合法或非法——彼等构造仅处于 DCG 的未定义域中，与本体系无涉。

DCG 不对欧氏几何、闵氏几何或其他公理体系内的构造作合法性评判。各体系相互独立，仅在自愿接受 DCG 公理约束时，DCG 方启动其合法性检验程序。DCG 无意替代任何既有几何学，亦不宣称自身为唯一真体系。

**体系定位**：DCG 是一个元几何框架，以五条操作性公理定义「合法几何闭包」的判定准则。满足公理的构造（欧氏、分形、微分几何、$q$-变形空间或高阶表示层）均被承认为合法几何对象；不满足者归于未定义域。

**操作顺序原则**：在 DCG 中，**构造优先于测量**。闭包必须首先通过公理 Ⅰ–Ⅳ 的检验成为合法几何对象，方可对其进行积分（公理 Ⅱ、Ⅴ 中出现的积分运算均在此原则下进行）。积分不参与构造，而是构造完成后的外边界性质提取。

**并置操作的积分规则**：
- 对于两个子域 $\mathcal{D}_1, \mathcal{D}_2$，其并置闭包 $\mathcal{D} = \mathcal{D}_1 \cup \mathcal{D}_2$ 的积分 $\mathcal{Z}_q(\mathcal{D})$ 仅在下述情形等于 $\mathcal{Z}_q(\mathcal{D}_1) + \mathcal{Z}_q(\mathcal{D}_2)$：
  - 当且仅当 $\mathcal{D}_1$ 与 $\mathcal{D}_2$ **严格分离**（$\mathcal{D}_1 \cap \mathcal{D}_2 = \varnothing$）。
- 对于**全等重叠、内部相交、边界贴合**三种情形，必须先合并、消去共享内面（或重叠区域），再在合并后的外边界上计算积分；不得先分块积分再相加作为并置闭包的积分。

---

### 1.2　公理 Ⅰ（紧致性）

设几何闭包 $\mathcal{O} = \langle \Gamma \rangle$ 为单一不可分构造。若满足
$$
\operatorname{conv}(\Gamma) = \mathcal{O}, \qquad \dim(\mathcal{O}) = \operatorname{rank}(\Gamma),
$$
且缩半操作 $H(\mathcal{O})$ 完全包含于 $\mathcal{O}$ 内部、不触及外部参照系，则称 $\mathcal{O}$ 是**紧致闭包**。

---

### 1.3　公理 Ⅱ（互易性）

设 $\mathcal{O} = \langle \Gamma \rangle$ 为紧致闭包。若其生成规则天然含有分支结构，则缩半操作 $H$ 自然诱导两个互补派生子闭包 $\mathcal{O}_{+}$、$\mathcal{O}_{-}$。若对规范内在度量族 $\{\mathcal{M}_k\}$ 中每一非平凡度量均有
$$
\mathcal{M}_k(\mathcal{O}_{+}) = \mathcal{M}_k(\mathcal{O}_{-}),
$$
则称 $\mathcal{O}$ 具有**自身互易性**。无明确分裂机制的闭包不声称互易，不通过本公理检验。

---

### 1.4　公理 Ⅲ（可逆性）

设递归生成序列 $\{\mathcal{O}_n\}_{n=0}^{\infty}$ 由映射 $\gamma_n \colon \mathcal{O}_n \to \mathcal{O}_{n+1}$ 给出。若对每个 $n$ 存在逆映射 $\gamma_n^{-1} \colon \mathcal{O}_{n+1} \to \mathcal{O}_n$ 满足
$$
\gamma_n^{-1} \circ \gamma_n = \mathrm{id}_{\mathcal{O}_n}, \qquad 
\gamma_n \circ \gamma_n^{-1} = \mathrm{id}_{\mathcal{O}_{n+1}},
$$
且 $\gamma_n^{-1}$ 保持生成基、度量族与配对结构，则称该递归是**可逆的**。极限闭包 $\mathcal{O}_\infty$ 称为**可逆无穷闭包**，若 $\gamma_\infty^{-1}$ 存在并保持结构。

---

### 1.5　公理 Ⅳ（一致性）

设 $\{\mathcal{O}_1, \dots, \mathcal{O}_m\}$ 为有限个紧致且互易的几何闭包。若其并置闭包 $\bigcup_{k=1}^{m} \mathcal{O}_k$ 在合并空间内仍紧致，则称该族是**一致的**。一致族的并置本身是一个合法几何闭包，可参与更高阶的一致族构造。

---

### 1.6　公理 Ⅴ（Dixmier 迹可积性与解析配对）

设 $\mathcal{O}_\infty$ 为可逆无穷闭包，其由**参数簇**驱动的非线性填充关联映射 $\Phi_{\text{nl}}$ 所诱导的局域配对泛函 $\omega(t)$ 不满足绝对可数可和性。若由递归结构自然诱导的尺度算子 $\mathcal{D}$ 使 $\omega(t)$ 属于紧算子理想 $\mathcal{L}^{(1,\infty)}$（Macaev 理想），则称该闭包具有 **Dixmier 可积性**。其全局积分算子定义为 Dixmier 奇异迹：
$$
\operatorname{Tr}_\omega(\mathcal{A}) := \lim_{N \to \infty} \frac{1}{\log N} \sum_{n=1}^N \lambda_n(\mathcal{A}),
$$
其中 $\lambda_n$ 为紧算子的特征值序列。该迹为公理 Ⅱ 中互易性判据在参数化扭曲情形下的唯一合法延拓。

**解析配对等价表述（梅林对偶）**：上述 Dixmier 迹等价于梅林生成函数在 $s=0$ 处的有限部分提取：
$$
\operatorname{Tr}_\omega = \text{F.P.}_{s=0}\, \Phi_{\mathcal{O}}(s), \qquad
\Phi_{\mathcal{O}}(s) := \sum_{n=1}^{\infty} \omega_n \cdot n^{-s},
$$
其中 $\text{F.P.}$ 表示洛朗展开的常数项系数。该等价性将公理 Ⅴ 转化为可操作的解析延拓判据。

> **总则**  
> 凡未通过上述五重公理检验的构造，DCG 不授予其几何闭包的合法性。

## 可选公理 Ⅵ & Ⅶ

**总则**：公理 Ⅵ、Ⅶ 均为可选约束。未声明者，DCG 仅依公理 Ⅰ–Ⅴ 检验；声明但核验失败者，仅撤销对应可选标签（「唯一」/「复合」），若通过公理Ⅰ–Ⅴ检验不剥夺基础合法地位。

---

### 可选公理 Ⅵ（唯一性）

**声明**：声明者须通过下述检验方可获授「唯一几何闭包」资格。

**定义 Ⅵ.1（自生递归序列）**  
设合法闭包 $\mathcal{O}$ 存在递归序列 $\{\mathcal{O}_n\}$ 收敛至 $\mathcal{O}$。若满足：

1. **生成内禀性**：$\gamma_n = \Gamma\big(\operatorname{Int}(\mathcal{O}_n)\big)$，$\Gamma$ 不依赖外部数据；
2. **参数闭合性**：$\operatorname{Int}(\mathcal{O}) = \lim\limits_{n\to\infty} \operatorname{Int}(\mathcal{O}_n)$，极限与序列选取无关（取 Kuratowski 集合极限）。

则称 $\{\mathcal{O}_n\}$ 为 $\mathcal{O}$ 的**自生递归序列**。

**公理 Ⅵ（唯一性）**  
若 $\mathcal{O}$ 存在自生递归序列，且任意两条满足 $\mathcal{O}_n = \mathcal{O}'_n$（$\forall n \in \mathbb{N}$），则称 $\mathcal{O}$ 满足唯一性。

**等价表述**：闭包是自身内生生成规则下唯一确定的极限几何体（不区分参数化形式）。

**实例**：科赫雪花通过；圆与椭圆闭包通常不声明；多组迭代规则收敛同一点集的变斜率填充体不通过。

---

### 可选公理 Ⅶ（复合性）

**声明**：凡涉及两个及以上合法闭包通过递归复合生成新闭包者，若声明本公理，须显式提供完整证明。

**定义 Ⅶ.1（复合递归生成）**  
设 $\{\mathcal{O}^{(i)}\}_{i=1}^m$（$m \ge 2$）已通过公理 Ⅰ–Ⅴ。若新闭包 $\mathcal{O}$ 由递归序列 $\{\mathcal{O}_n\}$ 生成，且每一 $\mathcal{O}_n$ 构造涉及至少两个 $\mathcal{O}^{(i)}$ 的并置、交叠、粘合或替换，$\gamma_n$ 显式依赖 $\operatorname{Int}(\mathcal{O}^{(i)})$，则称 $\mathcal{O}$ 为**复合递归生成**。

**公理 Ⅶ（复合性）**  
若 $\mathcal{O}$ 复合递归生成，须满足：

1. **组件合法性**：每个 $\mathcal{O}^{(i)}$ 独立通过公理 Ⅰ–Ⅴ；若声明 Ⅵ 则须同时满足。
2. **自证完备性**：显式证明 $\mathcal{O}$ 满足公理 Ⅰ–Ⅴ，溯源至组件结构，刻画复合操作对边界与内禀参数的改动。
3. **非平凡性**：复合合法性非组件合法性的简单推论，须独立重验。
4. **交叉约束（Ⅵ 联动）**：若同时声明 Ⅵ，还须证明：
   - (a) 任意自生递归序列对应的复合分解方式在等价意义下唯一；
   - (b) 不存在两组互不等价的组件分解方案生成同一极限闭包。

**实例**：柱状提升仅声明 Ⅶ 时完成张量分解证明即可，叠加 Ⅵ 则追加分解唯一性证明；混合曲率闭包仅声明 Ⅶ 时提供接续度规证明，叠加 Ⅵ 则论证无第二种拆分；两科赫雪花静态并置不触发 Ⅶ。

---

**声明组合与检验范围**：

| 声明组合 | 检验范围 |
|:---|:---|
| 未声明任一 | 仅公理 Ⅰ–Ⅴ |
| 仅 Ⅵ | 自生递归序列唯一性 |
| 仅 Ⅶ | 复合生成条件 1–3 |
| Ⅵ + Ⅶ | Ⅵ 全部 + Ⅶ 全部四条（含交叉约束） |

## 第二部分　几何本体论

### 2.1　全域（Universal Domain）

DCG 预设三条相互垂直的无限刻度线——横轴 $\mathcal{X}$、纵轴 $\mathcal{Y}$ 与竖轴 $\mathcal{Z}$，其刻度全体为实数集 $\mathbb{R}$。三轴张成背景地址空间
$$
\mathbb{R}^3 := \mathcal{X} \times \mathcal{Y} \times \mathcal{Z}.
$$

该空间**仅用于为几何点提供参照地址**，**不参与任何闭合构造**；其本身无界、无端点、无度量，不构成几何闭包。称 $\mathbb{R}^3$ 为 DCG 的**全域**。

全域固定不变。它的唯一功能是作为索引库，供后续从中切出有界子域。

---

### 2.2　子域（Sub-domain）与几何闭包

从全域 $\mathbb{R}^3$ 中取出三轴上的有界实数区间 $I_X, I_Y, I_Z \subset \mathbb{R}$，构成**约束域**
$$
\mathcal{D} := I_X \times I_Y \times I_Z \subset \mathbb{R}^3.
$$

若 $\mathcal{D}$ 满足缩半条件——即存在内部点 $c$ 使 $H_c(\mathcal{D}) \subsetneq \operatorname{int}(\mathcal{D})$——则称 $\mathcal{D}$ 为**合法子域**。由合法子域经复定坐标映射得到的三维内涵空间
$$
\boxed{\mathcal{S}}_q(\mathcal{D}) := \{\,\zeta_q(P) \mid P\in\mathcal{D}\,\}
$$
称为一个 **$q$-空间闭包**。单个闭包的承载流形是三维有限维的，不包含任何"无穷"。

---

### 2.3　域变量（Domain Variable）与空间域（Spatial Domain）

**域变量**：用于刻画合法子域 $\mathcal{D}$ 在 $\mathbb{R}^3$ 中"切法"的独立参数。

当子域取为轴对齐盒状时，域变量为六个区间端点：
$$
(a_X, b_X, a_Y, b_Y, a_Z, b_Z) \in \mathbb{R}^6,
$$
满足 $a_X < b_X$，$a_Y < b_Y$，$a_Z < b_Z$。更一般地，子域的边界可由可逆函数描述，此时域变量为函数空间中的元素。

**空间域**：设 $\mathfrak{D}$ 为 $\mathbb{R}^3$ 中全体合法子域的集合：
$$
\mathfrak{D} := \{\, \mathcal{D} \subset \mathbb{R}^3 \mid \mathcal{D} \text{ 为合法子域} \,\}.
$$
对每个 $\mathcal{D} \in \mathfrak{D}$，存在唯一的 $q$-空间闭包 $\boxed{\mathcal{S}}_q(\mathcal{D})$。全体闭包在商映射 $\bar{\zeta}_q$、广义轮换 $R_q$ 与外延积分算子 $\bar{\mathcal{Z}}_q$ 的全局一致性下粘结，构成一个整体：
$$
\boxed{\mathcal{S}}_q^{\,\text{dom}} := \bigcup_{\mathcal{D} \in \mathfrak{D}} \boxed{\mathcal{S}}_q(\mathcal{D}).
$$
称 $\boxed{\mathcal{S}}_q^{\,\text{dom}}$ 为 **$q$-空间域**。空间域的完整公理化见第六部分。

---

### 2.4　维度的严格区分

DCG 中"维"一词具有两重截然不同的含义：

| 维度类型 | 定义 | DCG 中的取值 |
|:---|:---|:---|
| **几何拓扑维数** | 承载闭包的空间流形的拓扑维数 | 恒为 3；商去 $\ker\zeta_q$ 后诱导的复结构实维数为 2 |
| **域变量维数** | 刻画子域"切法"的独立参数个数 | 由 $\mathfrak{D}$ 的选取决定；有限或无穷 |

**核心结论**：
- 几何闭包自身永远是有限维的（$\dim \leq 3$）。
- 空间域可以是无穷维的——当域变量的数目为无穷时。这种"无穷"仅代表**切法的选择无穷**，绝不代表底流形 $\mathbb{R}^3$ 的拓扑维数增加。

---

### 2.5 固有参数与派生参数

**定义 2.5.1（固有参数，Intrinsic Parameter）**

设 DCG 中已明确定义的实体层级为：

$$
\text{全域} \;\subset\; \text{空间域} \;\subset\; \text{子域} \;\subset\; \text{几何闭包} \;\subset\; \text{表示层}.
$$

对于任意实体 $E$ 属于上述层级之一，参数 $\xi$ 称为 $E$ 的**固有参数**，当且仅当 $\xi$ 满足以下条件之一：

1. $\xi$ 直接参与 $E$ 的构造定义，改变 $\xi$ 会导致 $E$ 的几何结构发生改变；
2. $\xi$ 是 $E$ 所属层级公理检验所需的内禀量。

固有参数的全体记作 $\mathrm{Int}(E)$。

**实例**：
- 子域的固有参数：区间端点 $(a_X, b_X, a_Y, b_Y, a_Z, b_Z)$，边界函数；
- 几何闭包的固有参数：曲率参数 $q$，填充斜率 $k_x, k_y$，递归深度 $n$；
- 表示层的固有参数：等差参数 $d$，伽罗瓦表示 $[\rho]$。

**定义 2.5.2（派生参数，Derived Parameter）**

对于实体 $E$，参数 $\tau$ 称为 $E$ 的**派生参数**，当且仅当：

1. $\tau$ 完全由 $E$ 的固有参数和 $E$ 上的构造操作（积分、度量、判定算法等）唯一确定；
2. $\tau$ 不满足固有参数的条件——改变 $\tau$ 不会改变 $E$ 的几何结构。

派生参数的全体记作 $\mathrm{Der}(E)$。

**实例**：
- 积分值 $\mathcal{Z}_q(\mathcal{O})$；
- 散度 $\operatorname{Div}_{\mathfrak{D}}(\mathcal{D})$、旋度 $\operatorname{Curl}_{\mathfrak{D}}(\mathcal{D})$、流度 $\operatorname{Flux}_{\mathfrak{D}}(\mathcal{D})$；
- 梅林生成函数 $\Phi_{\mathcal{O}}(s)$；
- 障碍类 $\mathrm{Ob}_{12}$；
- 互易性比值。

**公设 2.5.3（本体论层级隔离原则）**

设 $\tau \in \mathrm{Der}(E)$，其中 $E$ 属于层级 $\mathcal{L}$。则：

$$
\boxed{\tau \text{ 的语义解释域严格限定在 } \mathcal{L} \text{ 之内。}}
$$

即：

1. $\tau$ **不得**被解释为任何其他层级实体的固有参数；
2. $\tau$ **不得**反向参与 $E$ 的构造修改；
3. $\tau$ **不得**跨级用于推导其他层级实体的本体属性。

违反上述任意一条，均构成**范畴错误（category error）**，所得命题在 DCG 内无真值，归于未定义域。

**推论 2.5.4（时间参数与构造性无关禁令）**

设 $n$ 为几何闭包 $\mathcal{O}$ 的递归深度，满足 $n \in \mathrm{Int}(\mathcal{O})$。  
设 $\{ \mathcal{O}_n \}$ 为由 $n$ 索引的递归构造序列，$\{\mathcal{D}_i\}$ 为空间域中由域变量索引的子域族。

时间参数 $t$ 定义为外部观察者附加于序列 $\{ \mathcal{O}_n \}$ 或 $\{\mathcal{D}_i\}$ 的**描述性标签**，满足：

1. $t \in \mathrm{Der}(\mathcal{O})$；
2. $t$ 与 $n$ 之间**不存在任何构造性关联**：
   $$
   \boxed{\nexists\; f: \mathbb{N} \to \mathbb{R},\; f \in \mathrm{Aut}(\mathcal{O}),\; \text{使得 } t \equiv f(n).}
   $$
   $n$ 是 $\mathcal{O}$ 的内禀固有参数，$t$ 是观察者的事后赋值；
3. $t$ 的引入不参与 $\mathcal{O}$ 的公理检验，不改变 $\mathcal{O}$ 的几何结构；
4. 禁止由 $t$ 的存在推导 $\mathcal{O}$ 自身的流动性、膨胀性或收缩性：
   $$
   \text{流动性}(\mathcal{O}) \notin \{\text{真}, \text{假}\}, \quad \text{流动性}(\mathcal{O}) \in \text{未定义域}.
   $$

---

### 2.6　与标准泛函分析的对照

标准泛函分析中的“无穷维”来源于函数空间的像值自由度，底空间固定，拓扑维数不增加。DCG 的空间域与之形成精确对偶：空间中的“点”是几何点 $P \in \mathcal{D} \subset \mathbb{R}^3$，无穷仅来自子域的“切法”变化（域变量），底空间 $\mathbb{R}^3$ 固定不变。

## 第三部分　二维 DCG

### 3.1　操作基底

二维 DCG 是三维构造在 $z=0$ 平面上的限制。其背景为 $\mathcal{X}\mathcal{Y}$-平面，全域为 $\mathbb{R}^2$，合法子域为有界区间 $I_X \times I_Y$。

#### 3.1.1　直线递归闭包

给定两点 $P_1=(x_1,y_1)$、$P_2=(x_2,y_2)$（非重合），定义直线递归闭包为所有仿射组合：
$$
L(P_1, P_2) := \{\, (1-t)P_1 + tP_2 \mid t \in \mathbb{R} \,\}.
$$
此为无界闭包，仅当受约束限制时方有资格进入合法构造。

#### 3.1.2　非直线偏移

取定非常数偏移函数 $\delta(t)$（$t\in\mathbb{R}$）及方向 $\vec{w}\not\parallel (P_2-P_1)$，定义非直线递归生成元：
$$
\gamma_{\text{non}}(t) := P_1 + t(P_2-P_1) + \delta(t)\,\vec{w}.
$$
$\delta(t)$ 必须**取定**且支持回溯：存在确定的逆操作使路径从任意点还原至前驱生成元。不满足该条件的 $\delta$ 被视为非法。

#### 3.1.3　替换规则

生成映射 $\gamma \colon \mathcal{O}_n \to \mathcal{O}_{n+1}$ 将生成元替换为一组首尾相接的子路径，必须满足：
- **端点不变性**：替换路径的起点与终点与原生成元端点重合；
- **逆映射可定义**：从子路径可唯一识别父生成元；
- **度量按固定比例收缩**。

#### 3.1.4　旋转算子

旋转算子 $R_\theta$ 必须由体系内蕴构造给出。$R_{60^\circ}$ 由等边三角形第三顶点的几何构造定义：以两点 $A,B$ 为圆心、$|AB|$ 为半径作两圆，两圆交点即为第三顶点，由此确定 $R_{60^\circ}$ 的线性作用。

---

### 3.2　二维几何元素

#### 3.2.1　点（判定滤子）

设潜在交点池 $\mathfrak{B} := X\times Y \subset \mathbb{R}^2$。给定判定映射 $\Phi \colon \mathfrak{B} \to \mathbb{R}$，点集定义为零纤维
$$
\mathcal{L} := \Phi^{-1}(\{0\}) = \{\, (x,y)\in\mathfrak{B} \mid \Phi(x,y)=0 \,\}.
$$
仅当 $\Phi$ 为可测函数且 $\mathcal{L}$ 能纳入后续公理检验时，$\mathcal{L}$ 方被承认为几何点集。

#### 3.2.2　线

合法线闭包为直线递归闭包或非直线递归生成元经紧致性、互易性、可逆性筛选后的闭包。无界直线因不可逆而被排斥。

#### 3.2.3　面

设四条递归路径 $X_1, X_2, Y_1, Y_2$ 分别连接点 $A,B$、$C,D$、$A,C$、$B,D$，端点满足相交条件，形成首尾闭合环：
$$
Y_1 \xrightarrow{A} X_1 \xrightarrow{B} Y_2 \xrightarrow{D} X_2 \xrightarrow{C} Y_1.
$$
称此框架为一个**平面** $\mathbb{P} := \mathrm{Frame}(X_1,X_2,Y_1,Y_2)$。当某条路径退化为零长度时，三角形作为合法退化情形自动包含。

#### 3.2.4　圆（对数配对闭包）

设 $\Gamma = \{P_1, P_2, \dots, P_n\} \subset \mathbb{R}^2 \setminus \{(0,0)\}$ 为有限非零点集。固定参照点 $P_1 \in \Gamma$，定义**对数配对映射**：
$$
\Lambda_{P_1} \colon\ \Gamma \times \Gamma \to \mathbb{R}, \qquad (P_i, P_j) \mapsto \ln\|P_i + P_j\|.
$$
取定阈值 $\ell_0 \in \mathbb{R}$，构造平面点集：
$$
\mathcal{C}_{\ell_0} := \left\{\, Q \in \mathbb{R}^2 \setminus \{(0,0)\} \;\middle|\; \ln\|Q + P_1\| = \ell_0 \,\right\}.
$$
称 $\mathcal{C}_{\ell_0}$ 为由点集 $\Gamma$ 与参照点 $P_1$ 生成的**对数配对圆**。等价于以 $-P_1$ 为圆心、半径 $R=e^{\ell_0}$ 的欧氏圆。

该闭包由等式约束 $\Phi_{\ell_0}(Q) := \ln\|Q + P_1\| - \ell_0 = 0$ 确定，符合 §3.2.1 判定滤子之形式。其紧致、互易且可逆，通过全部合法性检验。

#### 3.2.5　平滑椭圆闭包

在 DCG 中，椭圆作为满足焦点约束与参数化光滑性的特殊闭包出现，是圆闭包向非常数焦距的自然推广，也是 $\mathcal{E}$-闭包分段协议的光滑极限。

**定义数据**：
- 焦点对 $F_1, F_2 \in \mathbb{R}^2$，$F_1 \neq F_2$。
- 主半轴长 $a \in \mathbb{R}^+$，满足 $a > c = \frac{1}{2}\|F_2 - F_1\|$。
- 方向向量 $\vec{w} \not\parallel (F_2 - F_1)$。
- 可选曲率参数 $q < 3/4$。

**隐式约束与参数化**：椭圆定义为点到两焦点距离之和为常数：
$$
\|Q - F_1\| + \|Q - F_2\| = 2a.
$$
令中心 $C = \frac{F_1+F_2}{2}$，主方向 $\vec{u} = \frac{F_2-F_1}{\|F_2-F_1\|}$，副方向 $\vec{v}$ 为 $\vec{u}$ 逆时针旋转 $90^\circ$，半短轴 $b = \sqrt{a^2 - c^2}$，则显式光滑参数化为
$$
\gamma(t) = C + a \cos(2\pi t)\,\vec{u} + b \sin(2\pi t)\,\vec{v}, \quad t \in [0,1].
$$
满足 $\gamma(0)=\gamma(1)$，$C^\infty$ 光滑。

**与圆及 $\mathcal{E}$-闭包的统一**：
- 当 $F_1 = F_2$ 时，$c=0$，退化为圆（取 $\ell_0 = \ln a$，$P_1 = -F_1$）。
- 作为 $\mathcal{E}$-闭包的光滑极限：若控制多边形为椭圆外切多边形，偏移函数 $\delta(t)$ 取为椭圆弧与弦的凸度差，则分段构造强收敛至上述 $\gamma(t)$。
- 曲率参数 $q$ 通过椭圆模 $\tau$ 与 $k^2 = 1 - b^2/a^2$ 关联，参数化可视为 Jacobi 椭圆函数的实截面，使椭圆在 $\zeta_q$ 下成为椭圆曲线的实嵌入。

**合法性**：紧致（有界闭集）、互易（上下半椭圆弧长相等）、可逆（$\gamma$ 为单射闭曲线）、一致，且不涉及参数簇发散，无需公理 Ⅴ，为完全合法闭包。

---

### 3.3　积分算子 $\mathcal{J}$ 与互易性判据

设 $(\gamma, T)$ 为闭包 $\mathcal{O}$ 的生成表示。定义局域配对泛函 $\langle\cdot,\cdot\rangle_{\mathcal{O}}$，局域贡献元 $\omega(t):=\langle\gamma(t), \mathrm{d}\gamma(t)\rangle_{\mathcal{O}}$。积分算子为沿生成元序列的代数累加：
$$
\mathcal{J}(\mathcal{O}) := \bigoplus_{t\in T} \omega(t),
$$
其中 $\bigoplus$ 表可数累加。该算子不依赖外部连续统测度，是几何内部的不变量。

**积分算子的边界性**：积分算子 $\mathcal{J}(\mathcal{O})$（以及其三维版本 $\mathcal{Z}_q$、$\bar{\mathcal{Z}}_q$）仅在合法闭包 $\mathcal{O}$ 的**外边界**上执行累加。内部结构（包括递归生成的历史、子闭包之间的贴合面、胞腔细分的内壁）不直接贡献于积分值。

**并置闭包的积分操作流程**：
1. 先判断两子域的拓扑关系（命题 6.1）。
2. 若为**严格分离**（$\mathcal{D}_1 \cap \mathcal{D}_2 = \varnothing$）：无需合并，直接分别积分再相加：
   $$
   \mathcal{Z}_q(\mathcal{D}_1 \cup \mathcal{D}_2) = \mathcal{Z}_q(\mathcal{D}_1) + \mathcal{Z}_q(\mathcal{D}_2).
   $$
3. 若为**全等重叠、内部相交、边界贴合**三者之一：必须**先合并**，消除共享内面（或重叠区域），使合并后的闭包通过公理 Ⅳ 检验，然后**仅在合并后的外边界上**计算积分：
   $$
   \mathcal{Z}_q(\mathcal{D}_1 \cup \mathcal{D}_2) \neq \mathcal{Z}_q(\mathcal{D}_1) + \mathcal{Z}_q(\mathcal{D}_2).
   $$

公理 Ⅱ 等价于积分等式：
$$
\mathcal{J}(\mathcal{O}_{+}) = \mathcal{J}(\mathcal{O}_{-}).
$$

若递归生成 $\gamma_n \colon \mathcal{O}_n \to \mathcal{O}_{n+1}$ 可逆（公理 Ⅲ），且所有有限闭包均互易（公理 Ⅱ），则极限闭包 $\mathcal{O}_\infty$ 必互易。

---

### 3.4　递归无穷闭包的合法性

递归无穷闭包 $\mathcal{O}_\infty = \lim_{n\to\infty}\mathcal{O}_n$ 的合法性由以下三则联合保证：
1. **归纳法**：基底 $\mathcal{O}_0$ 合法，且合法性在每一步生成中被保持。
2. **可逆性**（公理 Ⅲ）：每一步存在结构保持逆映射，确保任意有限步可回溯至基底。
3. **康托尔对角线法**：极限点被包含于初始紧致闭包内，不逃逸。

无界边界因不可逆，无法满足上述条件，故无穷无界闭包不合法。


## 第四部分　三维 DCG 与 $q$-变形空间

### 4.1　三维基底与内蕴曲率参数 $q$

三维 DCG 的几何本体承载于全域 $\mathbb{R}^3$ 的合法子域上。三轴地位等价，任意轮换 $\mathcal{X}\to\mathcal{Y}\to\mathcal{Z}\to\mathcal{X}$ 保持基底结构不变。

设 $q \in \mathbb{R}$ 为空间闭包的内蕴曲率参数：
- $q=0$：平直欧氏空间；
- $q>0$：椭圆型常曲率空间；
- $q<0$：双曲型（闵氏）常曲率空间。

---

### 4.2　$q$-变形代数元 $\omega_q$ 与复定坐标

定义 $\omega_q$ 为方程的根：
$$
\omega_q^2 + \omega_q + (1-q) = 0, \quad \omega_q \neq 1.
$$
显式解：
$$
\omega_q = \frac{-1 \pm \sqrt{4q-3}}{2}.
$$
$q=0$ 时退化为标准三次单位根 $\omega$。

对任意点 $P=(x,y,z)\in\mathbb{R}^3$，定义**内涵复定坐标**：
$$
\zeta_q(P) := x + y\,\omega_q + z\,\omega_q^2 \in \mathbb{R}(\omega_q).
$$

线性映射 $\zeta_q: \mathbb{R}^3 \to \mathbb{R}(\omega_q)$ 的核 $\ker\zeta_q$ 是一个一维子空间。商去该核，得到自然同构：
$$
\bar{\zeta}_q: \mathbb{R}^3 / \ker\zeta_q \xrightarrow{\;\sim\;} \mathbb{R}(\omega_q) \cong \mathbb{C} \quad (q<3/4).
$$
此商空间的实维数为 2，它作为复平面 $\mathbb{C}$ 的代数表示层出现，编码了三维循环对称性，但绝不引入第四个独立的空间坐标轴。称 $\bar{\zeta}_q$ 为**外延复定不变量**。

---

### 4.3　广义轮换算子 $R_q$

正确的 $q$-依赖轮换算子 $R_q$ 定义为满足下式的唯一线性映射：
$$
\zeta_q\bigl(R_q(P)\bigr) = \omega_q \cdot \zeta_q(P).
$$

展开 $\zeta_q$，解得：
$$
\boxed{R_q(x,y,z) = \bigl(z(1-q),\;\; x + zq,\;\; y\bigr)}.
$$

**性质**：

- $R_q$ 是 $\mathbb{R}^3$ 上的线性同胚，保持 $\ker\zeta_q$ 不变，并在商空间 $\mathbb{R}^3 / \ker\zeta_q$ 上诱导乘以 $\omega_q$ 的旋转变换。
- $R_q$ 的三次迭代在商空间上对应于乘以 $\omega_q^3$ 的变换，其中 $\omega_q^3 = (1-q) + q\,\omega_q$。当且仅当 $q = 0$ 时，$R_q^3 = \mathrm{id}$ 在 $\mathbb{R}^3$ 上整体成立。
- 当 $q=0$ 时，$R_0(x,y,z) = (z,x,y)$，退化为标准循环置换。
- 对于任意 $P$，有 $\zeta_q(R_q(P)) = \omega_q\,\zeta_q(P)$。
- 雅可比矩阵 $J_{R_q}$ 的行列式为 $\det(J_{R_q}) = 1-q$。在完整 $\mathbb{R}^3$ 上不恒为 1；但在商空间 $\mathbb{R}^3/\ker\zeta_q$ 上，$R_q$ 诱导的线性变换保持体积（行列式为 1）。

---

### 4.4　雅可比-可逆性前置判据

$q$-变形坐标映射 $\zeta_q:\mathbb{R}^3\to\mathbb{R}(\omega_q)$ 必须保持**实秩为 2**。经计算：
$$
\operatorname{rank}_{\mathbb{R}}\{1,\omega_q,\omega_q^2\} =
\begin{cases}
2, & q < 3/4,\\
1, & q \geq 3/4.
\end{cases}
$$

因此**合法曲率参数区间为 $q < 3/4$**。$q=3/4$ 为临界退化点，$q>3/4$ 全部非法。此判据优先于所有后续公理检验。

---

### 4.5　三维几何闭包

#### 4.5.1　$q$-空间闭包

取合法子域 $\mathcal{D}$，经复定坐标映射得到的三维内涵空间
$$
\boxed{\mathcal{S}}_q(\mathcal{D}) := \{\,\zeta_q(P) \mid P\in\mathcal{D}\,\}
$$
为一个 **$q$-空间闭包**，记作 $\boxed{\mathcal{S}}_q = \langle \mathcal{D}, q \rangle$，要求 $q<3/4$。

#### 4.5.2　柱状提升

设 $\mathcal{O}^{(2)}$ 为一合法的二维几何闭包。选定 $\mathcal{Z}$ 轴有界区间 $I_Z$，定义**三维几何闭包**为代数积：
$$
\boxed{\mathcal{O}}^{(3)} := \mathcal{O}^{(2)} \times I_Z = \{\, (x,y,z) \mid (x,y)\in\mathcal{O}^{(2)},\, z\in I_Z \,\}.
$$

#### 4.5.3　填充拉伸

设 $\boxed{\mathcal{O}}^{(3)} = \mathcal{O}^{(2)} \times I_Z$ 为一柱状闭包。对于每个 $z\in I_Z$，指定一对**填充边界函数** $\alpha_z(x), \beta_z(x)$，满足 $\alpha_z(x)\le\beta_z(x)$ 且均为可逆映射。则**填充拉伸闭包** $\boxed{\mathcal{F}}$ 由下式定义：
$$
\boxed{\mathcal{F}} := \{\, (x,y,z)\in\boxed{\mathcal{S}}_q \mid \exists (x_0,y_0)\in\mathcal{O}^{(2)},\, z\in I_Z,\; y\in[\alpha_z(x_0),\beta_z(x_0)] \,\}.
$$
该操作满足可逆性：通过剥离填充层可唯一恢复至原始柱状闭包。

#### 4.5.4　三维闭包簇

有限个 $q$-空间闭包与三维几何闭包若在并置后整体仍满足紧致性，则称为一个**三维闭包簇** $\boxed{\mathfrak{C}}$。


## 第五部分　互易性与判据

### 5.1　复定互易性

#### 5.1.1　互补子闭包的三维构造

对于柱状提升闭包，若二维闭包诱导互补分裂 $\mathcal{O}^{(2)} = \mathcal{O}^{(2)}_{+} \cup \mathcal{O}^{(2)}_{-}$，则三维闭包自然继承此分裂：
$$
\boxed{\mathcal{O}}^{(3)}_{+} := \mathcal{O}^{(2)}_{+} \times I_Z, \qquad 
\boxed{\mathcal{O}}^{(3)}_{-} := \mathcal{O}^{(2)}_{-} \times I_Z.
$$

#### 5.1.2　内涵复定积分算子

对闭包 $\boxed{\mathcal{F}}$ 的生成路径 $\gamma$，定义**内涵复定积分算子**：
$$
\mathcal{Z}_q(\boxed{\mathcal{F}}) := \bigoplus_{t\in T} \langle\, \zeta_q(\gamma(t)),\; \mathrm{d}\zeta_q(\gamma(t)) \,\rangle_{\mathbb{R}(\omega_q)}.
$$

#### 5.1.3　复定互易性判据

设紧致三维闭包 $\boxed{\mathcal{F}}$ 的生成规则自然诱导互补子闭包 $\boxed{\mathcal{F}}_{+}$、$\boxed{\mathcal{F}}_{-}$。若
$$
\mathcal{Z}_q(\boxed{\mathcal{F}}_{+}) = \mathcal{Z}_q(\boxed{\mathcal{F}}_{-}),
$$
则称 $\boxed{\mathcal{F}}$ 具有**复定互易性**。

#### 5.1.4　柱状提升定理

若二维闭包 $\mathcal{O}^{(2)}$ 自身互易，且填充边界函数在互补区域上一致，则乘积闭包或填充拉伸闭包自动满足复定互易性。此时三维复定积分满足分解律：
$$
\mathcal{Z}_q\big(\boxed{\mathcal{F}}_{\pm}\big) 
\;=\; 
\mathcal{J}_q\big(\mathcal{O}^{(2)}_{\pm}\big) \;\otimes\; \Lambda(I_Z).
$$

**推论**：二维互易性沿纤维方向自然"粘贴"为三维互易性，且与曲率参数 $q<3/4$ 的具体取值无关。

---

### 5.2　线性斜率填充体与外延判据

设 $\Phi$ 引入线性交叉耦合项：
$$
\Phi(x,y,z) = \big(x + k_x z,\ y + k_y z,\ z\big),
$$
其中 $k_x,k_y \in \mathbb{R}$ 为**填充斜率参数**。称由该映射生成的实体 $\boxed{\mathcal{F}}_{k_x,k_y}$ 为**斜率填充体**。

展开 $\zeta_q$：
$$
\zeta_q(\Phi) = x + y\omega_q + z \cdot \Xi_q(k_x,k_y),
$$
其中
$$
\Xi_q(k_x,k_y) = (k_x + q - 1) + (k_y - 1)\omega_q.
$$

**外延复定判据**：对于斜率填充体，互易性判定改写为外延积分形式：
$$
\bar{\mathcal{Z}}_q(\boxed{\mathcal{F}}_{+}) = \bar{\mathcal{Z}}_q(\boxed{\mathcal{F}}_{-}),
$$
其中 $\bar{\mathcal{Z}}_q$ 为基于 $\bar{\zeta}_q$ 的外延复定积分算子：
$$
\bar{\mathcal{Z}}_q(\boxed{\mathcal{F}}) := \bigoplus_{t\in T} \langle\, \bar{\zeta}_q(\gamma(t)),\; \mathrm{d}\bar{\zeta}_q(\gamma(t)) \,\rangle_{\mathbb{R}(\omega_q)}.
$$

该判据自动吸收全部线性交叉耦合项，与斜率参数 $k_x,k_y$ 的取值无关。

---

### 5.3　参数簇

**定义（参数簇）**：在固定合法子域 $\mathcal{D}$ 内，若填充斜率参数依赖于递归深度或地址位置的函数 $k_x(z), k_y(z)$，则这些参数的全体构成一个**参数簇**。其参数空间可为无穷维函数空间，但每个具体取值仍作用于有限维底空间 $\mathbb{R}^3$ 上。

此时填充关联映射为：
$$
\Phi_{\text{nl}}(x,y,z) = \big(x + k_x(z) \cdot z,\ y + k_y(z) \cdot z,\ z\big).
$$

普通可数累加积分算子 $\bigoplus$ 在极限闭包 $\mathcal{O}_\infty$ 处可能发散。此时须启用公理 Ⅴ。

---

### 5.4　雅可比-梅林联合判据

设第 $n$ 层递归的雅可比矩阵为 $J_n$，定义扰动范数：
$$
\delta_n := \|J_n - I\|_{\text{op}}.
$$

| 扰动行为 | 梅林生成函数 | 判据 | 合法性 |
|:---|:---|:---|:---|
| $\sum \delta_n < \infty$（绝对收敛） | $\Phi(s)$ 在 $s=0$ 全纯 | 普通累加 $\bigoplus$ | ✅ 合法 |
| $\delta_n \sim 1/n$（对数发散） | $\Phi(s)=\zeta(s+1)$，$s=0$ 简单极点 | Dixmier 迹 / 梅林有限部分 | 🔄 公理 Ⅴ 启用 |
| $\delta_n \sim n^{-\alpha}$，$\alpha<1$ | $\Phi(s)$ 在 $s=0$ 非可去奇性 | Dixmier 迹失效 | ❌ 非法 |

**判定树**：$\alpha>1$ → 合法；$\alpha=1$ → 个案（公理 Ⅴ）；$\alpha<1$ → 非法。

**梅林乘法配对判据**：
$$
\langle\!\langle \mathcal{O}_+, \mathcal{O}_- \rangle\!\rangle_{\mathrm{Mellin}} := \text{F.P.}_{s=0} \left[ \Phi_{\mathcal{O}_+}(s) \cdot \overline{\Phi_{\mathcal{O}_-}(s)} \right].
$$
若此有限部分在互补分裂下对称，则基于参数簇的闭包通过互易性检验。


## 第六部分　空间域整体变换

### 6.1　概述

空间域整体变换是 DCG 中固定底空间 $\mathbb{R}^3$、仅改变子域"切法规则"的完整框架。它不增加坐标轴数量，不改变几何点的三个自由度，而是通过伽罗瓦表示层对子域族进行分类，并通过解陪域处理不同表示层之间的粘合障碍。

**核心原则**：底流形永远固定为 $\mathbb{R}^3$（拓扑维数 = 3）。无穷性完全来自子域切法规则的多样性及其相互之间的对称关系。

---

### 6.2　无穷谱系

DCG 中的「无穷」分布于四个独立层级。凡涉及无穷的构造，必须首先明确其所属层级。

| 层级 | 名称 | 来源 | 承载空间 | 控制手段 |
|:---|:---|:---|:---|:---|
| 一 | **递归无穷** | 迭代深度 $n \to \infty$ | 有限维 | 归纳法 + 可逆性回溯 |
| 二 | **参数无穷** | 函数空间自由度 | 有限维 | 公理 Ⅴ（Dixmier / 梅林） |
| 三 | **并置无穷** | 闭包数量 $\to \infty$ | 局部有限 | 紧致 exhaustion + 局部有限覆盖 |
| 四 | **表示无穷** | 伽罗瓦表示层数量 $\to \infty$ | 纤维丛结构 | 循环上同调 + 解陪域筛选 |

---

### 6.3　空间域的可逆相交拓扑

**定义 6.1（子域可逆邻域）**  
设 $\mathcal{D}_0 \in \mathfrak{D}$。对 $\mathcal{D} \in \mathfrak{D}$，若存在自同构
$$
\lambda_{\mathcal{D}} \in \operatorname{Aut}(\mathfrak{D})
$$
满足：

1. $\lambda_{\mathcal{D}}$ 保持缩半条件、保持 $\bar{\zeta}_q$ 代数相容性；
2. $\lambda_{\mathcal{D}}(\mathcal{D}_0) \cap \mathcal{D} \neq \varnothing$；
3. $\lambda_{\mathcal{D}}$ 在交集上限制为恒等映射：
   $$
   \lambda_{\mathcal{D}}|_{\mathcal{D}_0 \cap \mathcal{D}} = \mathrm{id},
   $$

则称 $\mathcal{D}$ 属于 $\mathcal{D}_0$ 的**可逆邻域**。

全体可逆邻域生成 $\mathfrak{D}$ 上的**可逆相交拓扑**。  
在此拓扑下，**子域族 $\mathfrak{D}$ 本身成为拓扑空间**。

> **关键区分**  
> 该拓扑是**以子域为点的集合拓扑**，不是 $\mathbb{R}^3$ 的欧氏子空间拓扑。

---

**命题 6.1（子域四型拓扑关系）**  
任意两个合法子域 $\mathcal{D}_1, \mathcal{D}_2 \in \mathfrak{D}$ 的拓扑关系被四类开集穷尽：

| 类型 | 条件 | 对应开集 |
|:---|:---|:---|
| 全等重叠 | $\mathcal{D}_1 = \mathcal{D}_2$ | 平凡单点邻域 |
| 内部相交 | $\operatorname{int}(\mathcal{D}_1) \cap \operatorname{int}(\mathcal{D}_2) \neq \varnothing$，边界无交 | 内部交集非空 |
| 边界贴合 | $\partial\mathcal{D}_1 \cap \partial\mathcal{D}_2 \neq \varnothing$ | 边界交集非空 |
| 严格分离 | $\mathcal{D}_1 \cap \mathcal{D}_2 = \varnothing$ | 补集邻域 |

四类关系互斥穷尽，构成空间域拓扑基。

---

### 6.4　空间域无穷维雅可比与域变量高斯判定

**定义 6.2（空间域雅可比 $J_{\mathrm{dom}}$）**  
设 $\{\mathcal{D}_i\}_{i \in I} \subset \mathfrak{D}$ 为一族合法子域，$\{\xi_j\}$ 为**域变量**（边界模态、区间端点、外形自由度）。

定义空间域雅可比矩阵
$$
(J_{\mathrm{dom}})_{ij}
:=
\left.\frac{\partial \zeta_q(\mathcal{D}_i)}{\partial \xi_j}\right|_{\xi = \xi_0}.
$$

> **严格层级隔离**  
> $J_{\mathrm{dom}}$ **仅作用于域变量（层级三：换框自由度）**。  
> **不包含** $k_x(z), k_y(z)$ 内部形变参数（层级二：单框内自由度），  
> 亦区别于递归细分雅可比 $J_n$（层级一：迭代扰动）。

空间域雅可比的行列式与整体奇性由 **公理 Ⅴ Dixmier 迹 / 梅林有限部分** 规范，无需绝对收敛。

---

**算法 6.1（空间域高斯消元判定算法）**

**输入**：$\mathcal{D}_1, \mathcal{D}_2 \in \mathfrak{D}$，对应行向量 $r_1, r_2$。

**步骤 1**：构造无穷维增广矩阵
$$
M = \begin{pmatrix} r_1 \\ r_2 \end{pmatrix}.
$$

**步骤 2**：可逆行化简（**有限支撑约束**）  
所有行变换限于**有限支撑扰动**：即每一行操作（行交换、行倍加、行缩放）仅涉及有限个非零域变量列。  
此限制确保行变换在可逆相交拓扑下为同胚，保持邻域结构的局部性，不会将远处无关切法混入局部相交判定。

> **良定义性注记**  
> 可逆相交拓扑的开集由有限次基本初等变换生成。有限支撑约束恰好对应这一拓扑基的有限交性质，保证消元过程中每一中间步骤都停留在同一可逆邻域内，不会跨域跳变。因此算法输出与消元路径选择无关，良定义。

**步骤 3**：主元分离  
设主元列为 $p_1, p_2$：

- 若 $p_1 = p_2$ 且分量全同 → **重叠候选**；
- 若 $p_1 \neq p_2$ → 构造分离参数：
  $$
  \Lambda_{\text{sep}} =
  \begin{cases}
  \displaystyle \frac{(J_{\mathrm{dom}})_{2,p_2}}{(J_{\mathrm{dom}})_{1,p_1}} & \text{乘法型} \\[6pt]
  (J_{\mathrm{dom}})_{2,p_2} - (J_{\mathrm{dom}})_{1,p_1} & \text{加法型}
  \end{cases}
  $$

**步骤 4**：粗判定

- $\Lambda_{\text{sep}} \neq \mathbf{1}\,/\,0$ → **严格分离**；
- $\Lambda_{\text{sep}} = \mathbf{1}\,/\,0$ → 进入步骤 5。

**步骤 5**：边界/相交精细判别  
取边界诱导行向量（仅取有限支撑投影）：

- 边界行线性相关 → **边界贴合**；
- 边界行线性无关 → **纯内部相交**。

**输出**：$\{\text{重叠},\ \text{相交},\ \text{边界贴合},\ \text{分离}\}$ 唯一结果。

---

**命题 6.2（有限支撑行变换拓扑相容性）**  
设 $\mathcal{N}(\mathcal{D}_0)$ 为 $\mathcal{D}_0$ 在可逆相交拓扑下的任一开邻域。算法 6.1 中的任意行变换若限于有限支撑扰动，则它将 $\mathcal{N}(\mathcal{D}_0)$ 同胚地映到另一个可逆邻域 $\mathcal{N}(\mathcal{D}_0')$。

*证明概要*：有限支撑行变换仅作用有限个域变量，每个域变量对应一个可逆重参量化方向，有限个方向的连续形变仍落在同一可逆邻域的有限交中，逆变换同样有限支撑，故为同胚。

---

### 6.5　相交可逆无损还原定理

**定理 6.1（相交可逆性与无损拼接）**  
若 $\mathcal{D}_1, \mathcal{D}_2$ 满足相交条件，则存在全局可逆映射
$$
\gamma_{\cap}^{-1}
:
\zeta_q(\mathcal{D}_1 \cap \mathcal{D}_2)
\to
\zeta_q(\mathcal{D}_1) \times \zeta_q(\mathcal{D}_2)
$$
满足：

1. **互易结构守恒**：保持公理 Ⅱ 互补配对；
2. **算法可构造**：$\gamma_{\cap}^{-1}$ 完全由空间域高斯行变换显式生成；
3. **局部兼容**：对任意子域 $\mathcal{D}' \subset \mathcal{D}_1 \cap \mathcal{D}_2$，
   $$
   \gamma_{\cap}^{-1}|_{\mathcal{D}'}
   =
   \gamma_1^{-1}|_{\mathcal{D}'} \oplus \gamma_2^{-1}|_{\mathcal{D}'}.
   $$

*证明概要*：两子域各自具有公理 Ⅱ 互补分裂结构。相交条件保证交叉正负配对在雅可比化简后完全等价。交集商结构可唯一分解为两子域各自逆结构的直和，故信息无丢失、粘合无矛盾。

> **核心推论**  
> 全体子域闭包 $\{\boxed{\mathcal{S}}_q(\mathcal{D})\}$ **可以无矛盾代数粘结**，是 $\boxed{\mathcal{S}}_q^{\infty}$ 整体存在性的严格理论基础。

**拼接后积分的唯一性**：定理 6.1 保证了拼接后的闭包 $\mathcal{O}_1 \cup \mathcal{O}_2$ 在通过公理 Ⅳ 检验后是一个合法的几何闭包。其积分
$$
\mathcal{Z}_q(\mathcal{O}_1 \cup \mathcal{O}_2)
$$
仅在合并后的**外边界**上执行，共享面（内部面）已被消除。因此：

- 当 $\mathcal{O}_1$ 与 $\mathcal{O}_2$ **严格分离**时：
  $$
  \mathcal{Z}_q(\mathcal{O}_1 \cup \mathcal{O}_2) = \mathcal{Z}_q(\mathcal{O}_1) + \mathcal{Z}_q(\mathcal{O}_2).
  $$

- 在其他三种拓扑关系下：
  $$
  \mathcal{Z}_q(\mathcal{O}_1 \cup \mathcal{O}_2) \neq \mathcal{Z}_q(\mathcal{O}_1) + \mathcal{Z}_q(\mathcal{O}_2).
  $$

该性质不是 DCG 的缺陷，而是其构造逻辑的一贯性：积分是闭包的外在属性，构造完成后，边界唯一决定积分值。

---

### 6.6　空间域的平凡纤维丛结构

**命题 6.3（空间域整体平凡丛）**  
DCG 三维无穷空间域具有标准平凡纤维丛结构：
$$
\boxed{\mathcal{S}}_q^{\infty}
\;\cong\;
\mathfrak{D} \times \mathbb{R}(\omega_q).
$$

- **底空间**：$\mathfrak{D}$（带可逆相交拓扑的合法子域族）；
- **标准纤维**：$\mathbb{R}(\omega_q)$ 二维复 $q$-商结构；
- **每根纤维**：单一子域对应的完整 $q$-闭包 $\boxed{\mathcal{S}}_q(\mathcal{D})$。

*证明概要*：对任意 $\mathcal{D}$，存在可逆重参量化 $\rho_{\mathcal{D}} : \mathcal{D} \to \mathcal{D}_0$ 与商映射 $\zeta_q$ 可交换。所有局部平凡化通过可逆邻域拓扑连通、无绕、无全局障碍，因此整体平凡。

---

> **物理几何图景（最终标准表述）**  
> - **底流形永远固定 $\mathbb{R}^3$（拓扑维数 = 3）**；  
> - **无穷性完全来自底空间 $\mathfrak{D}$ 的域变量自由度**。  
> 
> 完美区分：
> - 传统泛函：固定定义域、跑函数（轴固定、态无穷）；
> - DCG 空间域：固定映射规则、跑定义域（态固定、框无穷）。

---

### 6.7　空间域的相对紧致化

**定义 6.3（相对紧致空间域）**  
给定半径界 $R > 0$，定义有界约束子域全体：
$$
\boxed{\mathcal{S}}_q^{\infty, R}
:=
\left\{
\mathcal{D} \in \mathfrak{D}
\;\middle|\;
\sup_{P \in \mathcal{D}} \|P\| \le R
\right\}.
$$

该集合在**可逆相交拓扑**下紧致。

全体 $\{\boxed{\mathcal{S}}_q^{\infty, R}\}_{R > 0}$ 构成 $\boxed{\mathcal{S}}_q^{\infty}$ 的**紧致穷尽递增族**。无穷空间域由有限紧致版本极限生成，全程满足 DCG 紧致化程序，无发散、无逃逸、结构稳定。

---

### 6.8　度量完备化

设 $\gamma_n$ 的收缩比为 $r_n \in (0,1)$，即
$$
\mathcal{M}_{k}(\mathcal{O}_{n+1}) = r_n \cdot \mathcal{M}_{k}(\mathcal{O}_n).
$$

**命题**：若 $\sum_{n=0}^{\infty} (1-r_n) < \infty$，则存在唯一的内蕴度量 $\mathcal{M}_\infty$ 定义于 $\mathcal{O}_\infty$ 上，使极限闭包在此度量下完备。

---

### 6.9　紧致化程序

DCG 不排斥非紧致空间，但要求任何涉及无穷的构造必须通过**紧致化程序**。

**单点紧致化闭包**：设 $\mathcal{O}$ 为局部紧致但非紧致的构造。若存在 $\infty \notin \mathcal{O}$ 及紧致邻域系 $\{K_n\}$，满足 exhaustion、紧致性、零贡献条件，则称 $\mathcal{O}_\bullet := \mathcal{O} \cup \{\infty\}$ 为 $\mathcal{O}$ 的 DCG 单点紧致化闭包。

**紧致化与互易性**：若 $\mathcal{O}$ 自身具有互补分裂，则紧致化继承分裂，且由于 $\omega(\infty)=0$，互易性在紧致化下保持。

---

### 6.10　公理 Ⅴ 扩展判定算法（含解陪域恢复）

**算法 6.2（公理 Ⅴ 扩展判定算法）**

**输入**：参数簇，配对泛函序列 $\{\omega_n\}$，以及可能涉及的伽罗瓦表示层 $[\rho_1], [\rho_2]$。

**步骤 1**：计算扰动范数 $\delta_n := \|J_n - I\|_{\text{op}}$。

**步骤 2**：判断 $\sum \delta_n$ 收敛性。
- 若收敛 → **合法**（公理 Ⅰ–Ⅳ）。
- 若发散 → 进入步骤 3。

**步骤 3**：构造梅林生成函数 $\Phi_+(s), \Phi_-(s)$。

**步骤 4**：计算 $\Delta\Phi(s) := \Phi_+(s) - \Phi_-(s)$，提取 $s=0$ 处洛朗展开，得到障碍类 $\mathrm{Ob}$。

**步骤 5**：判定奇性：
- 全纯 → **合法**；
- 简单极点留数为零 → **条件合法**；
- 非可去奇性 → **进入步骤 6**。

**步骤 6（解陪域恢复检查）**：
- 设障碍类为 $\mathrm{Ob}_{12} \neq 0$（对应于表示层 $[\rho_1], [\rho_2]$ 之间的粘合相位不匹配）。
- 检查是否存在解陪域 $\mathfrak{R}_{12}$，使得提升后的障碍 $\mathrm{Ob}_{12}^{(\eta)} = 0$。
- 检验 $\mathfrak{R}_{12}$ 是否满足公理 Ⅰ–Ⅳ：
  - 若 $\mathfrak{R}_{12}$ 合法 → **合法（经由解陪域恢复）**；
  - 若不存在这样的 $\mathfrak{R}_{12}$，或构造出的 $\mathfrak{R}_{12}$ 不通过公理 Ⅰ–Ⅳ → **非法（不可补救）**。

**输出**：$\{\text{合法},\ \text{合法（经解陪域）},\ \text{非法}\}$。

---

**元定理（解陪域完备性）**：  
对于 DCG 中任意非零循环上同调障碍 $\mathrm{Ob}_{12}$，若存在一个紧致、可逆、互易、一致的解陪域 $\mathfrak{R}_{12}$ 使得 $\mathrm{Ob}_{12}^{(\eta)} = 0$，则该障碍不构成合法性失效，而仅构成 **"构造复杂度增加"** 的标志。

---

### 6.11　表示层筛选器与解陪域构造

**定义 6.4（伽罗瓦表示层）**  
设 $G = \operatorname{Gal}(K/\mathbb{Q})$ 为 $q$-变形代数元定义域的伽罗瓦群。对每个二维实表示 $[\rho] \in \operatorname{Rep}_2(G)$，定义子域族
$$
\mathfrak{D}_{[\rho]} := \{\, \mathcal{D} \in \mathfrak{D} \mid \mathcal{D} \text{ 的切法规则兼容表示 } [\rho] \,\}.
$$

不同表示层可以有不同的缩半规则、不同的紧致化程序，以及不同的域变量参数化方式。

**定义 6.5（循环上同调粘合障碍）**  
设 $[\rho_1], [\rho_2]$ 为两个表示层，其纤维在交集 $\mathcal{D} \in \mathfrak{D}_{[\rho_1]} \cap \mathfrak{D}_{[\rho_2]}$ 上的过渡函数为 $g_{12}$。定义循环上同调障碍类：
$$
\mathrm{Ob}_{12} := [g_{12}] \in HC^1(\mathcal{D}, \operatorname{Iso}(F_1, F_2)).
$$

- 若 $\mathrm{Ob}_{12} = 0$：两层可无碍粘合；
- 若 $\mathrm{Ob}_{12} \neq 0$：需引入解陪域。

**定义 6.6（解陪域）**  
设 $[\rho_1], [\rho_2] \in \operatorname{Rep}_2(G)$，其过渡函数 $g_{12}$ 的循环上同调障碍为 $\mathrm{Ob}_{12} \neq 0$。称纤维层 $\mathfrak{R}_{12}$ 为 $g_{12}$ 的**解陪域**，若存在相位松弛参数 $\eta \in \mathfrak{R}_{12}$ 使得提升过渡函数
$$
\tilde{g}_{12}^{(\eta)} := \eta \circ g_{12}
$$
满足
$$
\mathrm{Ob}_{12}^{(\eta)} = 0.
$$

此时复合构造
$$
[\rho_1] \xrightarrow{\tilde{g}_{12}^{(\eta)}} \mathfrak{R}_{12} \xrightarrow{\tilde{g}_{21}^{(\eta)}} [\rho_2]
$$
为合法 DCG 闭包。

---

**命题 6.4（解陪域存在性）**  
若 $g_{12}$ 的障碍类 $[\mathrm{Ob}_{12}] \in HC^1(A)$ 为**挠元（torsion）**，则存在**有限维解陪域** $\mathfrak{R}_{12}$。  
若 $[\mathrm{Ob}_{12}]$ 非挠，则解陪域维数无穷，须启用**公理 Ⅴ** 的 Dixmier 迹判据控制其合法性。

---

**命题 6.5（解陪域的 DCG 合法性）**  
解陪域 $\mathfrak{R}_{12}$ 本身必须满足公理 Ⅰ–Ⅳ 的全部要求。若 $\mathfrak{R}_{12}$ 为无穷维，则其迹由公理 Ⅴ 规范。

---

**核心推论**：  
DCG 的合法粘合判定分为三种情形：

| 情形 | 障碍类 | 解陪域 | 合法性 |
|:---|:---|:---|:---|
| 零障碍 | $[\mathrm{Ob}_{12}] = 0$ | 不需要 | 合法 |
| 挠障碍 | $[\mathrm{Ob}_{12}]$ 为挠元 | 有限维，显式可构造 | 合法 |
| 非挠障碍 | $[\mathrm{Ob}_{12}]$ 非挠 | 无穷维，须 Dixmier 迹收敛 | 条件合法（由公理 Ⅴ 判定） |
| 不可解析 | $\mathrm{Ob}_{12}$ 不可嵌入任何合法解陪域 | 不存在 | 非法 |

此判据将「障碍非零 → 非法」扩展为「障碍非零 → 检查是否存在解陪域 → 若存在则合法」。

---

**复合纤维丛的完整结构**：
$$
\boxed{\mathcal{S}}_q^{(\infty,\text{Gal},\text{Res})}
=
\bigcup_{[\rho] \in \operatorname{Rep}_2(G)}
\left(
\mathfrak{D}_{[\rho]}
\times
\mathbb{R}(\omega_q^{(\rho)})
\times
\bigotimes_{\langle i,j \rangle} \mathfrak{R}_{ij}
\right)
\Big/ \sim_{\text{cyclic}},
$$
其中粘合 $\sim_{\text{cyclic}}$ 在所有非零障碍被对应解陪域吸收后总是合法的。


## 第七部分　无穷维度的统一构造

### 7.1　概述

本节处理 DCG 中坐标轴数量 $n > 3$ 的推广。它与第六部分「空间域整体变换」有本质区别：

- **空间域整体变换**：固定 $\mathbb{R}^3$，改变子域切法规则（域变量的无穷）。
- **无穷维度构造**：增加坐标轴数量 $n$（几何拓扑维数的拓展）。

两者共享同一套表示层筛选器（伽罗瓦群 + 循环上同调 + 解陪域），但作用于不同的层次。

---

### 7.2　核心构造：等比等差双模递推

**定义 7.1（$n$ 轴代数元系）**  
设 $n \geq 3$，取 $n$ 次单位根 $r = e^{2\pi i / n}$，等差参数 $d \in \mathbb{R}$，定义 $n$ 个代数元
$$
\boxed{w_k = r^k + k \cdot d, \qquad k = 0, 1, \dots, n-1.}
$$

---

**性质 7.1（轮换闭合）**  
$w_k$ 在 $n$ 次轮换后整体平移 $n \cdot d$。轮换闭合条件应理解为**模 $n d$ 平移等价类**：
$$
w_{k+n} = w_k + n \cdot d \;\equiv\; w_k \pmod{n d}.
$$

解陪域 $\mathfrak{R}_{12}$ 中的 $\eta$ 参数正是用来在不同表示层之间调节这一平移等价类的匹配。

---

### 7.3　广义轮换算子 $R_q^{(n)}$

**定义 7.2（$n$ 轴轮换算子）**  
设 $P = (x_0, x_1, \dots, x_{n-1}) \in \mathbb{R}^n$，定义
$$
\zeta_q^{(n)}(P) := \sum_{k=0}^{n-1} x_k \, w_k.
$$

广义轮换算子 $R_q^{(n)}$ 是 $\mathbb{R}^n$ 上的唯一线性映射，满足
$$
\zeta_q^{(n)}\bigl(R_q^{(n)}(P)\bigr) = w_1 \cdot \zeta_q^{(n)}(P).
$$

**性质 7.3（$n$ 重循环对称骨架）**  
设 $r = e^{2\pi i / n}$，则 $\{r^k\}_{k=0}^{n-1}$ 在复平面上构成一个 $n$ 重循环对称骨架：$r^{k+1}$ 由 $r^k$ 旋转 $2\pi/n$ 得到，$(R_q^{(n)})^n = \mathrm{id}$。这一循环结构是 $w_k$ 的普遍代数基础，直接保证 $n$ 轴轮换在任意 $n \geq 3$ 时满足公理 Ⅱ 的互易性——互补配对由旋转对称性自然诱导，与等差平移项 $k \cdot d$ 无关。

当 $n$ 恰好为 $3$ 的倍数时，上述 $n$ 重骨架在局部三元组 $\{r^k, r^{k+1}, r^{k+2}\}$ 上额外呈现 $120^\circ$ 相位差。此 $120^\circ$ 结构是 $n$ 重循环对称在 $n=3m$ 时的内禀特例，并非公理 Ⅱ 的普遍必要条件。无论 $n$ 是否为 $3$ 的倍数，$n$ 重循环骨架本身已足以支撑互易性判定。

$w_k = r^k + k \cdot d$ 的任意相邻三元组是在此 $n$ 重循环骨架上叠加了线性平移项 $k \cdot d$ 的结果。平移项改变的是三元组在复平面上的绝对位置，不改变骨架的循环对称性。因此公理 Ⅱ 的判定只检验 $r^k$ 部分，平移部分由解陪域 $\mathfrak{R}_{12}$ 的 $\eta$ 参数独立调节。

---

### 7.4　实秩判据与合法曲率参数

**定理 7.1（$n$ 轴实秩判据）**  
对于 $w_k = r^k + k \cdot d$，其中 $r = e^{2\pi i / n}$，$d \in \mathbb{R}$：
$$
\operatorname{rank}_{\mathbb{R}}\{w_0, w_1, \dots, w_{n-1}\} =
\begin{cases}
2, & d \neq 0 \text{ 或 } n \text{ 为 } 3 \text{ 的倍数时 } r \text{ 为三次单位根},\\
1, & \text{退化情形}.
\end{cases}
$$

**推论**：对任意 $n \geq 3$，存在 $d \neq 0$ 使得 $\{w_k\}$ 实秩为 2，即 $n$ 轴 $q$-变形空间合法嵌入复平面，保持 $q < 3/4$ 判据。

---

### 7.5　解陪域的显式形式

**定理 7.2（解陪域的显示构造）**  
设 $[\rho_1], [\rho_2]$ 对应的等差参数分别为 $d_1, d_2$。其解陪域为
$$
\mathfrak{R}_{12} :=
\begin{cases}
\{\, d_1 + m \cdot \Delta \mid m = 0, 1, \dots, M-1 \,\}, & \text{若 } [\mathrm{Ob}_{12}] \text{ 为挠元，阶数为 } M, \\
\mathbb{R}, & \text{若 } [\mathrm{Ob}_{12}] \text{ 非挠},
\end{cases}
$$
其中 $\Delta$ 是基本周期。

提升过渡函数为
$$
\tilde{g}_{12}^{(\eta)} := \eta \circ g_{12}, \qquad \eta \in \mathfrak{R}_{12},
$$
即用 $\eta$ 替换 $d_1$ 和 $d_2$，使得两层在粘合处共享相同的等差平移参数。

*证明概要*：$\eta$ 的调整仅改变等差平移项 $k \cdot \eta$，不触及 $r^k$ 的 $n$ 重循环对称骨架（后者是公理 Ⅱ 的保障来源），因此 $\mathrm{Ob}_{12}^{(\eta)} = 0$ 对所有 $\eta \in \mathfrak{R}_{12}$ 成立。

---

### 7.6　统一构造的完整公式

$$
\boxed{
\begin{aligned}
w_k &= r^k + k \cdot d, \quad r = e^{2\pi i / n}, \quad d \in \mathbb{R}, \\
R_q^{(n)} &: w_k \mapsto w_{k+1}, \quad \text{在 } r^k \text{ 上诱导纯旋转；} k \cdot d \text{ 部分随轮换平移，其不匹配由 } \mathfrak{R}_{12} \text{ 吸收}, \\
\mathfrak{R}_{12} &= \{\, \eta \in \mathbb{R} \mid \eta \text{ 替换 } d \text{ 的平移参数} \,\}, \\
\boxed{\mathcal{S}}_q^{(n)} &= \left\{ \sum_{k=0}^{n-1} x_k w_k \;\middle|\; (x_0, \dots, x_{n-1}) \in \mathcal{D} \subset \mathbb{R}^n \right\}, \quad q < 3/4.
\end{aligned}
}
$$


## 第八部分　矢量几何与复合构造

### 8.1　矢量场

设几何闭包 $\mathcal{O}$ 的生成路径为 $\gamma(t)$，切向量 $\mathbf{T}(t) := \frac{\mathrm{d}\gamma}{\mathrm{d}t}$，法向量 $\mathbf{n}(t)$ 由 $\mathbf{T}(t)$ 旋转 $90^\circ$ 得到。定义**内蕴矢量场**：
$$
\mathbf{V}(\gamma(t)) := v(t)\,\mathbf{T}(t) + w(t)\,\mathbf{n}(t).
$$

**离散散度与旋度**：
$$
\operatorname{curl}(\mathcal{O},\mathbf{V}) := \bigoplus_{t\in T} \langle \mathbf{V}(\gamma(t)), \mathbf{T}(t) \rangle_{\mathcal{O}},
$$
$$
\operatorname{div}(\mathcal{O},\mathbf{V}) := \bigoplus_{t\in T} \langle \mathbf{V}(\gamma(t)), \mathbf{n}(t) \rangle_{\mathcal{O}}.
$$

---

### 8.2　矢量互易性

若闭包具有互补分裂，且矢量场在分裂下保持对称，则称满足**矢量互易性**：
$$
\operatorname{curl}(\mathcal{O}_+,\mathbf{V}) = \operatorname{curl}(\mathcal{O}_-,\mathbf{V}), \qquad
\operatorname{div}(\mathcal{O}_+,\mathbf{V}) = \operatorname{div}(\mathcal{O}_-,\mathbf{V}).
$$

---

### 8.3　混合曲率空间

设 $\{\boxed{\mathcal{S}}_{q_i}\}_{i=1}^m$ 为有限个 $q$-空间闭包，$q_i<3/4$ 但不必相等。若并置后仍紧致，则称**混合闭包簇**。相邻空间交界面必须配备**接续度规** $\mathfrak{g}_{ij}$，满足保可逆性、度量连续、Dixmier 一致性。

---

### 8.4　分形曲面填充体

在柱状提升的每一层 $z\in I_Z$ 上，用二维闭包自身的递归分形边界作为填充边界函数：
$$
\alpha_z(x) := f_z(x), \qquad \beta_z(x) := g_z(x),
$$
称所得闭包为**分形曲面填充体**。其公理检验由基底闭包的合法性继承。


## 第九部分　物理场论嵌入

### 9.1　离散场量与作用量

定义**离散场量**：
$$
\varphi \colon \Gamma \to \mathbb{R} \text{ 或 } \mathbb{R}(\omega_q).
$$

离散作用量：
$$
S(\varphi; \mathcal{O}) := \mathcal{Z}_q(\mathcal{O}) \cdot \bigoplus_{t\in T} \langle \varphi(\gamma(t)), \Delta \varphi(\gamma(t)) \rangle_{\mathcal{O}},
$$
其中 $\Delta$ 为离散拉普拉斯算子。

---

### 9.2　运动方程与场论互易性

取变分 $\delta S = 0$ 给出**离散运动方程**。若 $S(\varphi; \mathcal{O}_+) = S(\varphi; \mathcal{O}_-)$，则称具有**场论互易性**。


## 附录

### 附录 A　约束的必要性

若放弃公理 Ⅳ，每构造一个复合几何体均需手工验证端点重合、边界无溢出、性质传递，体系不封闭。五重公理将验证前置至个体闭包资质审查，使 DCG 成为自足的规范几何系统。

---

### 附录 B　$q$ 参数合法性速查表

| $q$ 范围 | $\omega_q$ 类型 | 实秩 | $R_q$ | 合法性 |
|:---|:---|:---|:---|:---|
| $q < 3/4$ | 复数 | 2 | $R_q(x,y,z) = (z(1-q),\, x+zq,\, y)$ | ✅ 合法 |
| $q = 3/4$ | 实数（$-1/2$） | 1 | 退化 | ❌ 临界退化，非法 |
| $q > 3/4$ | 实数 | 1 | 退化 | ❌ 非法 |

---

### 附录 C　层级拓展对照表

| 扭曲类型 | 填充关联特征 | 积分/判据算子 | 互易性判据 | 适用公理 |
|:---|:---|:---|:---|:---|
| 标准柱状 | $k_x=k_y=0$ | 内涵 $\mathcal{Z}_q$ | $\mathcal{Z}_q(\mathcal{F}_+)=\mathcal{Z}_q(\mathcal{F}_-)$ | Ⅰ–Ⅳ |
| 线性斜率 | $k_x,k_y$ 常数 | 外延 $\bar{\mathcal{Z}}_q$ | $\bar{\mathcal{Z}}_q(\mathcal{F}_+)=\bar{\mathcal{Z}}_q(\mathcal{F}_-)$ | Ⅰ–Ⅳ |
| 参数簇 | $k_x(z), k_y(z)$ 变化 | Dixmier / 梅林 | $\text{F.P.}_{s=0}[\Phi_+ - \Phi_-]=0$ | Ⅰ–Ⅴ |
| 混合曲率 | 不同 $q$ 拼接 | 接续度规 $\mathfrak{g}_{ij}$ | 界面度量连续 | Ⅰ–Ⅴ |
| 分形曲面 | 分形边界沿 $Z$ 递归 | $\mathcal{Z}_q$ 或 Dixmier | 基底互易 + 插值一致 | Ⅰ–Ⅴ |
| 场论嵌入 | 离散作用量 | $\mathcal{Z}_q$ 或 Dixmier | $S(\mathcal{O}_+)=S(\mathcal{O}_-)$ | Ⅰ–Ⅴ |

---

### 附录 D　术语表

| 术语 | 定义 |
|:---|:---|
| 全域（Universal Domain） | 背景地址空间 $\mathbb{R}^3$，不参与构造 |
| 子域（Sub-domain） | 从全域切出的有界合法约束域 $\mathcal{D}$ |
| 域变量（Domain Variable） | 刻画子域"切法"的独立参数 |
| 空间域（Spatial Domain） | 全体合法子域闭包的并置，记作 $\boxed{\mathcal{S}}_q^{\,\text{dom}}$ |
| 空间域整体变换 | 固定 $\mathbb{R}^3$，通过伽罗瓦表示层与解陪域改变子域切法规则的完整框架 |
| 几何闭包（Geometric Closure） | 通过五重公理检验的合法构造 |
| 参数簇（Parameter Family） | 依赖连续参数的函数族，参数空间可为无穷维 |
| 解陪域（Resolution Fiber Domain） | 吸收循环上同调粘合障碍的中间纤维层，由等差参数 $d$ 的调节空间构成 |
| 等比等差构造 | $n$ 轴情形下 $w_k = r^k + k \cdot d$ 的统一代数结构 |
| 有限支撑扰动 | 行变换仅涉及有限个非零域变量列，保证拓扑同胚与良定义性 |

---

### 附录 E　构造实例（概要）

#### E.1　科赫雪花

初始线段 $P_1=(0,0), P_2=(1,0)$，生成映射 $\gamma_n$ 替换为四段：
$$
C = A + \tfrac{1}{3}(B-A),\; E = A + \tfrac{2}{3}(B-A),\; D = C + R_{60^\circ}(E-C).
$$
收缩比 $1/3$，积分（面积）收敛至 $A_\infty = \frac{8}{5}\cdot\frac{\sqrt{3}}{4}$，互补分裂面积恒等。满足公理 Ⅰ–Ⅳ。

#### E.2　克莱因瓶

初始正方形 $[0,1]^2$，$2\times2$ 胞腔细分并赋予粘合 $(0,y)\sim(1,y)$，$(x,0)\sim(1-x,1)$。单次细分完成，互易性由对角线分割保证，面积各半。合法有限闭包。

#### E.3　三维柱状提升

对二维互易闭包 $\mathcal{K}^{(2)}$，取 $I_Z=[0,H]$，闭包 $\boxed{\mathcal{K}}^{(3)} = \mathcal{K}^{(2)} \times I_Z$。积分分解：
$$
\mathcal{Z}_q(\boxed{\mathcal{K}}^{(3)}_{\pm}) = \mathcal{J}_q(\mathcal{K}^{(2)}_{\pm}) \otimes \Lambda(I_Z),
$$
互易自动继承。

#### E.4　常数斜率填充体

$\Phi(x,y,z) = (x + k_x z,\; y + k_y z,\; z)$，逆映射显式存在。外延积分吸收 $\Xi_q$，互易性与 $k_x,k_y$ 无关。通过 Ⅰ–Ⅳ。

#### E.5　变斜率填充体

$k_x(z), k_y(z)$ 为光滑函数，如 $k_x(z)=0.3\sin(2\pi z)$。扰动范数 $\delta_n \sim 1/n$，梅林生成函数 $\Phi_{\pm}(s)$ 在 $s=0$ 有简单极点但留数为零，经 Dixmier 迹提取有限部分得 $\Delta\mathcal{Z}_q = 0$，互易成立。若斜率非光滑致 $\delta_n \sim n^{-\alpha},\alpha<1$ 则非法。适用公理 Ⅴ。

#### E.6　极度尖锐科赫雪花变体（二维）

初始线段，替换规则改为极度尖锐等腰三角形（顶角 \(\varepsilon \to 0\)，高/底比无界增长），收缩比仍为 \(1/3\)。该构造满足公理 Ⅰ–Ⅳ，但扰动范数 \(\delta_n \sim n^{-\alpha}\)（\(0<\alpha<1\)），级数幂次发散且非对数型，梅林生成函数在 \(s=0\) 处非可去奇性，Dixmier 迹失效，公理 Ⅴ 不通过；同时递归抹平可逆性因跨层尖刺缠绕失效，自生递归序列不存在，公理 Ⅵ 亦无法通过。二维单层判定为非法。

#### E.7　三维分层复合（合法性救赎）

将上述尖刺沿第三维分层安置，各层由 \(z\) 坐标隔离，避免跨层缠绕；引入 \(z\) 方向权重 \(w(z)\) 压制高层扰动，使有效扰动 \(\delta_n^{3D} = \delta_n \cdot w(z_n)\) 的累积收敛或对数发散，从而通过公理 Ⅴ。该三维复合构造须经公理 Ⅶ 核验：组件（各分层单元）满足 Ⅰ–Ⅳ，且须证明整体合法不能由组件合法直接推出（本案例中组件单层原本非法，更无继承问题）。若同时声明公理 Ⅵ，还需论证复合分解方式唯一性。

---

### 附录 F　与标准数值工具的对照

**核心区别**：
- 规则区域：DCG 与标准数值积分一致；
- 分形边界：DCG 原生支持递归顶点序列，标准工具无法直接处理；
- 边界贴合：DCG 通过公理 Ⅳ 合并消重，共享内面只计一次；标准工具需手动定向消除重复。

**边界贴合验证实例**：$\mathcal{D}_1$ 与 $\mathcal{D}_4$ 共享 $x=2$ 面，合并后该面成为内部面被消除，并置闭包的积分值不等于分块积分之和，严格符合 DCG 积分规则。

---

### 附录 G　通用符号

| 符号 | 含义 |
|:---|:---|
| $\bigoplus$ | 可数代数累加（不依赖测度） |
| $\text{F.P.}_{s=0}$ | 洛朗展开在 $s=0$ 处的常数项 |
| $\boxed{\mathcal{S}}_q$ | 带有曲率参数 $q$ 的空间闭包 |
| $\mathfrak{D}$ | 全体合法子域的集合 |

---

### 附录 H　版本更新记录

| 版本 | 日期 | 更新内容 |
|:---|:---|:---|
| DCG 1.0 | — | 四重公理体系；二维/三维闭包构造；$q$-变形空间；矢量几何初步框架。 |
| DCG 1.5 | — | 新增公理 Ⅴ；外延复定结构；斜率填充体；分层判定框架。 |
| DCG 1.6 | — | 梅林对偶强化；混合曲率接续度规；分形曲面填充体；物理场论嵌入；$q<3/4$ 前置判据。 |
| DCG 1.7 | — | 广义轮换算子 $R_q$ 修正；$\ker\zeta_q$ 商空间修正；雅可比-梅林判据整合；$\Xi_q$ 闭式公式。 |
| DCG 2.0 | — | 本体论重构；全域/子域/域变量/空间域四层概念；无穷闭包理论三级谱系；度量完备化；紧致化程序；公理 Ⅴ 解析判定算法。 |
| DCG 2.2 | — | 空间域完整公理化移至附录展望；无穷维度（$n>3$）移入展望；全局术语统一；大纲重组。 |
| DCG 2.5 | — | 新增第六部分「空间域整体变换」；无穷谱系扩展为四层级；公理 Ⅴ 扩展判定算法（含解陪域恢复）；伽罗瓦表示层筛选器；循环上同调粘合障碍。 |
| **DCG 3.0** | — | 新增第七部分「无穷维度的统一构造」，采用等比等差双模递推 $w_k = r^k + k \cdot d$ 给出任意 $n \geq 3$ 轴情形的完整代数构造；解陪域显式形式由等差参数 $d$ 自然给出；无穷维度从"展望"升格为"已完成构造"；修正 4.3 节 $R_q^3 = \mathrm{id}$ 的表述；修正 7.3 节性质 7.3，明确 120° 循环结构属于 $r^k$ 骨架而非 $w_k$。 |
| **DCG 3.1** | — | 在 1.1 声明中完善"操作顺序原则"，明确严格分离与其它三种情形的积分规则差异；在 3.3 节积分算子边界性中增加"并置闭包的积分操作流程"表格；在 6.5 节定理 6.1 后补充"拼接后积分的唯一性"精确表述；修正 4.3 节雅可比行列式表述，明确完整 $\mathbb{R}^3$ 上 $\det(J_{R_q}) = 1-q$ 而在商空间上保体积；附录 F 新增积分对比结论表，严格区分四种拓扑关系下的积分等式成立条件；全文修订编号及交叉引用。进一步修正 7.3 节性质 7.3，将 120° 局部对称明确降级为 $n=3m$ 时的内禀特例，确立 $n$ 重循环对称骨架为公理 Ⅱ 的普遍保障来源；新增 §3.2.5 平滑椭圆闭包。 |
| **DCG 3.2** | — | 新增可选公理 Ⅵ（唯一性）与可选公理 Ⅶ（复合性），二者均为可选约束，声明且核验失败者仅撤销对应标签，不影响公理 Ⅰ–Ⅴ 下的基础合法地位；公理 Ⅵ 定义自生递归序列，以生成内禀性与参数闭合性为判据；公理 Ⅶ 定义复合递归生成，要求组件合法性、自证完备性、非平凡性，并与公理 Ⅵ 设交叉约束（联动条款）；声明组合与检验范围明确为四种状态；实例更新涵盖柱状提升、混合曲率、科赫雪花并置等场景。 |
| **DCG 3.3** | — | 新增附录 E.6「极度尖锐科赫雪花变体（二维）」边界案例分析：该构造满足 Ⅰ–Ⅳ，但扰动范数呈幂次发散（$\alpha<1$），公理 Ⅴ 失效；同时递归抹平可逆性因跨层尖刺缠绕失效，公理 Ⅵ 亦不通过，二维单层判定非法。新增附录 E.7「三维分层复合——合法性救赎」：将尖刺沿第三维分层安置并引入 $z$ 方向权重 $w(z)$，有效扰动收敛或对数发散，通过公理 Ⅴ；触发公理 Ⅶ 复合性核验，组件满足 Ⅰ–Ⅳ，非平凡性要求自动满足；若同时声明公理 Ⅵ 则附加论证分解唯一性。此二案例展示了公理 Ⅴ 的判别力、公理 Ⅲ 的“参数可逆”与“构造可逆”之区分，以及公理 Ⅶ 在三维分层复合中的介入逻辑。 |