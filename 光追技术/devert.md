# DyVeRT 技术设计文档

**版本 1.4**

---

## 1. 摘要

DyVeRT 是一个面向纯动态场景的 1 SPP Monte Carlo 路径追踪渲染器。系统采用双层 BVH 空间加速结构，每帧全量重建顶层聚合；渲染阶段执行全采样逐像素路径追踪，不采用探针、泛化或任何降采样策略；后处理在单帧空域内完成，以跨域双边滤波（Cross-Bilateral Filter）为唯一降噪模块，配合 Reinhard 色调映射与自适应 Bloom。场景昼夜差异全部收进统一参数配置表，主循环零硬编码。

---

## 2. 系统架构

```
输入场景 + 人物模型
       │
       ▼
 自动相机系统（局部封闭法）
       │
       ▼
 双层 BVH 构建 / 动态更新
       │
       ▼
 全采样逐像素路径追踪（1 SPP，多光源，硬阴影，指数雾，雨雪）
       │
       ▼
 HDR 辐射度 + Albedo + Normal 缓冲区
       │
       ▼
 跨域双边滤波（空间 / 颜色 / 法线三权重联合）
       │
       ▼
 Reinhard 色调映射 + sRGB 校正
       │
       ▼
 多级高斯 Bloom（昼夜参数分离）
       │
       ▼
 最终 LDR 图像
```

---

## 3. 自动相机系统（局部封闭法）

### 3.1 设计目标

给定场景所有静态几何的包围盒顶点、球心、树中心以及动态人物中心，自动确定相机位姿，保证主体完整可见且避免贴脸穿模。

### 3.2 通用局部封闭相机公式

#### 步骤 1：空间分块

将顶点集 $V$ 均匀格子化，格宽 $h$ 取场景平均边长的 $2\sim3$ 倍：

$$
\mathbf{idx}(\mathbf{v}) = \left\lfloor \frac{\mathbf{v} - \mathbf{v}_{\min}}{h} \right\rfloor
$$

#### 步骤 2：局部团提取

对每个非空格子 $\mathbf{g}$，提取其 $3\times3\times3$ 邻域内的所有顶点：

$$
V_{\mathbf{g}} = \bigcup_{\|\mathbf{g}'-\mathbf{g}\|_\infty \le 1} \{\mathbf{v}\in V \mid \mathbf{idx}(\mathbf{v}) = \mathbf{g}'\}
$$

#### 步骤 3：局部封闭性检验（体积感判定）

计算该团三轴跨度：

$$
\Delta_j(\mathbf{g}) = \max_{\mathbf{v}\in V_{\mathbf{g}}} v_j - \min_{\mathbf{v}\in V_{\mathbf{g}}} v_j, \quad j\in\{x,y,z\}
$$

若 $\min(\Delta_x, \Delta_y, \Delta_z) > \epsilon$（$\epsilon \approx 0.05h$），则认为该团具备体积感，通过封闭检验。

#### 步骤 4：局部形心

$$
\mathbf{c}_{\mathbf{g}} = \frac{1}{|V_{\mathbf{g}}|}\sum_{\mathbf{v}\in V_{\mathbf{g}}}\mathbf{v}
$$

#### 步骤 5：选择最佳局部团

选择体量最大且通过封闭检验的团：

$$
\mathbf{g}^* = \arg\max_{\mathbf{g}} |V_{\mathbf{g}}| \cdot \mathbb{I}\!\left[\min_j \Delta_j(\mathbf{g}) > \epsilon\right]
$$

#### 步骤 6：相机外推

计算全局均值：

$$
\mathbf{o} = \frac{1}{|V|}\sum_{\mathbf{v}\in V}\mathbf{v}
$$

定义逃离方向：

$$
\hat{\mathbf{d}} = \frac{\mathbf{c}_{\mathbf{g}^*} - \mathbf{o}}{\|\mathbf{c}_{\mathbf{g}^*} - \mathbf{o}\|}
$$

相机位置沿该方向外推：

$$
\mathbf{c}_{\text{cam}} = \mathbf{c}_{\mathbf{g}^*} + \lambda \cdot \|\mathbf{v}_{\max} - \mathbf{v}_{\min}\| \cdot \hat{\mathbf{d}}, \quad \lambda \in [1.5,\; 2.5]
$$

视线方向始终指向局部团形心 $\mathbf{c}_{\mathbf{g}^*}$。

#### 兜底策略

若所有团均未通过体积检验（纯散点、线框、开放曲面）：

$$
\mathbf{c}_{\text{cam}} = \mathbf{o} + 2\|\mathbf{v}_{\max} - \mathbf{v}_{\min}\| \cdot \hat{\mathbf{z}}
$$

直接全局均值上方俯视。

### 3.3 简化实现（备注）

若依赖受限或场景规模较小，可采用等价简化：将全局包围盒中心视为局部团中心，用固定方向外推。

$$
\begin{aligned}
\mathbf{o} &= \frac{\mathbf{v}_{\min} + \mathbf{v}_{\max}}{2} \\
d_{\text{len}} &= \|\mathbf{v}_{\max} - \mathbf{v}_{\min}\| \\
\mathbf{c}_{\text{cam}} &= \mathbf{o} + 0.6 \cdot d_{\text{len}} \cdot \frac{\mathbf{d}}{\|\mathbf{d}\|}
\end{aligned}
$$

其中 $\mathbf{d} = (2.0,\ -3.0,\ 1.5)$，视线 $\mathbf{lookat} = \mathbf{o} + (0,0,0.5)$。

---

## 4. 核心数据结构

### 4.1 统一顶点池 $G$

$$
G \in \mathbb{R}^{N \times K}
$$

每行 $G_i$ 为顶点完整局部空间属性向量：

| 列 | 属性 | 量纲 | 说明 |
|---|---|---|---|
| 0‑2 | $(x,y,z)$ | [L] | 局部坐标 |
| 3‑5 | $(nx,ny,nz)$ | [1] | 局部法线（单位向量） |
| 6‑7 | $(u,v)$ | [1] | 纹理坐标 |
| 8‑10 | $(tx,ty,tz)$ | [1] | 切线（可选） |
| 11‑14 | $(r,g,b,a)$ | [1] | 顶点颜色 / Albedo |
| 15+ | $\text{mat\_id}$ 等 | [1] | 材质实例 ID、流形图表 ID |

### 4.2 变换表 $T$

$$
T = \{M_1, M_2, \dots, M_m\}, \quad M_j \in \mathbb{R}^{4 \times 4}
$$

每个对象对应一个 4×4 齐次矩阵。法线和切线用逆转置矩阵保持正交性：

$$
\mathbf{n}_{\text{world}} = \text{normalize}\left( (M_j^{-1})^{\mathsf{T}} \cdot \begin{bmatrix} nx \\ ny \\ nz \\ 0 \end{bmatrix} \right)
$$

### 4.3 全局索引缓冲 $I$

$$
I \in \mathbb{N}^{N_I}
$$

所有三角形索引统一为 32 位无符号整数。每个 Primitive 在 $I$ 中占据连续段。

### 4.4 几何对象句柄

```cpp
struct GeometryHandle {
    uint32_t vertex_start, vertex_count;
    uint32_t index_start,  index_count;
    uint32_t transform_id;
    AABB     local_aabb;
    uint32_t material_id;
    uint32_t primitive_id;
};
```

---

## 5. 双层 BVH 空间加速结构

### 5.1 局部 BVH（对象级）

构建于对象局部坐标系，叶节点存储三角形索引三元组 $(i_0, i_1, i_2)$，指向全局索引缓冲 $I$。

#### 5.1.1 自适应宽度（法线一致性驱动）

对节点 $N$ 包含的面片集 $\mathcal{P}_N$：

$$
\bar{\mathbf{n}}(N) = \frac{1}{|\mathcal{P}_N|} \sum_{i \in \mathcal{P}_N} \mathbf{n}_i
$$

$$
\sigma(N) = \max_{i \in \mathcal{P}_N} \arccos(\mathbf{n}_i \cdot \bar{\mathbf{n}}(N))
$$

子节点个数 $k$：

$$
k(N) = \text{clamp}\left( \text{round}\left( k_{\max} \cdot \left(1 - \frac{\sigma(N)}{\sigma_{\max}}\right)^{\gamma} \right),\ k_{\min},\ k_{\max} \right)
$$

| 参数 | 推荐值 |
|---|---|
| $k_{\min}$ | 2 |
| $k_{\max}$ | 8 |
| $\sigma_{\max}$ | $\pi/2$ |
| $\gamma$ | 1.5 |
| 叶子容量 | 6 |

#### 5.1.2 构建算法（LBVH + 自适应宽度分割）

1. 统计法线，计算 $\sigma$ 与 $k$。
2. Morton 编码排序面片中心。
3. 将已排序序列分割为 $k$ 个连续块。
4. 递归构建子树。

### 5.2 顶层 TLAS（全局聚合）

TLAS 叶节点只存轻量元数据：

$$
\mathcal{L}_j = \left(\text{AABB}_{\text{world}}^{(j)},\ \text{ptr}(\mathcal{B}_j),\ j\right)
$$

世界包围盒由局部根 AABB 经 $M_j$ 变换 8 个角点得到。

#### 5.2.1 每帧全量重建

纯动态下 refit 会产生历史残留。全量 LBVH 构建保证 TLAS 拓扑严格反映当前帧空间分布。复杂度 $O(M \log M)$，$M \le 10^4$ 时 $< 0.1\ \text{ms}$。

#### 5.2.2 可见性粗筛

TLAS 构建前用视锥体 6 平面与世界 AABB 做点积测试，提前剔除不可见对象。

---

## 6. 动态更新策略

### 6.1 形变对象（顶点位置变化）

1. 更新 $G$ 中该对象顶点的局部坐标及法线。
2. 重新计算三角形面法线。
3. 重建该对象的局部 BVH。
4. 重新计算世界 AABB。
5. 标记 TLAS 需要重建（或直接全量重建）。

复杂度 $O(n \log n)$，仅限该对象。

### 6.2 刚性变换（仅 $M_j$ 变化）

1. 计算新的世界 AABB（从局部根 AABB 变换）。
2. 更新 TLAS 对应叶节点的 AABB。
3. 局部 BVH 拓扑完全不动。

零几何重建开销。

---

## 7. 渲染管线（主路径）

### 7.1 射线生成（针孔相机）

给定相机参数，像素 $(i,j)$ 的射线：

$$
\mathbf{o} = \mathbf{e}
$$

$$
\mathbf{d} = \text{normalize}\left( \mathbf{w} + \frac{2i - W}{H} \cdot \tan(fov/2) \cdot \mathbf{r} + \frac{2j - H}{H} \cdot \tan(fov/2) \cdot \mathbf{u} \right)
$$

### 7.2 双层遍历

遍历 TLAS 获取候选对象；将射线变换到对象局部空间后遍历局部 BVH；世界距离由两点世界坐标计算，禁止非均匀缩放下直接使用局部 $t$。

### 7.3 水密三角形求交

采用 Woop 等人算法，投影至主轴平面，2D 精确叉积测试，$\epsilon = 2^{-20}$。

### 7.4 材质与光照

**漫反射（Lambertian）**：

$$
L_{\text{diffuse}} = \text{albedo} \cdot (1 - \text{metal}) \cdot \frac{\max(0,\ \mathbf{N}\cdot\mathbf{L})}{\pi}
$$

**镜面反射（Beckmann 法线分布近似）**：

$$
\begin{aligned}
\mathbf{H} &= \frac{\mathbf{L} + \mathbf{V}}{\|\mathbf{L} + \mathbf{V}\|} \\
p &= \max\!\left(\frac{2}{\text{rough}^2 + 10^{-3}} - 2,\; 1\right) \\
L_{\text{specular}} &= (1 - \text{rough}) \cdot \text{metal} \cdot (\max(0,\ \mathbf{N}\cdot\mathbf{H}))^{p}
\end{aligned}
$$

**光源衰减**：

$$
\text{att} = \frac{\text{intensity}}{\|\mathbf{L}\|^2 + 1}
$$

**硬阴影**：从着色点向光源方向检测遮挡，若被遮挡则漫反射项乘以 $0.03$。

**环境光**：

$$
L_{\text{ambient}} = \text{albedo} \cdot 0.02
$$

**噪声注入**：零均值高斯噪声，标准差 $\sigma \sqrt{|\text{color}| + 0.05}$，其中 $\sigma$ 为场景基准。

**指数雾**：

$$
C_{\text{final}} = C \cdot e^{-d \cdot t} + C_{\text{fog}} \cdot (1 - e^{-d \cdot t})
$$

**天气粒子**：雨雪在夜晚开启，直接叠加条纹或高斯光斑到 HDR 图像。

---

## 8. G-Buffer 输出

- **HDR 颜色**：线性空间辐射度。
- **Albedo**：无光照表面色。
- **Normal**：世界空间法向量（已归一化）。

---

## 9. 后处理管线

### 9.1 跨域双边滤波（Cross-Bilateral Filter）

对每个颜色通道独立执行。像素 $(x,y)$ 的邻域权重由三部分乘积构成：

$$
w(i,j) = \exp\!\left(-\frac{(i-x)^2+(j-y)^2}{2\sigma_s^2}\right)
\cdot \exp\!\left(-\frac{\|\mathbf{A}_{ij} - \mathbf{A}_{xy}\|^2}{2\sigma_r^2}\right)
\cdot \exp\!\left(-\frac{\|\mathbf{N}_{ij} - \mathbf{N}_{xy}\|^2}{2\sigma_n^2}\right)
$$

滤波结果：

$$
I'(x,y) = \frac{\sum_{(i,j)\in\Omega} w(i,j) \, I(i,j)}{\sum w(i,j)}
$$

其中 $\Omega$ 为以 $(x,y)$ 为中心的 $(2r+1)\times(2r+1)$ 窗口。

| 参数 | 符号 | 值 | 说明 |
|---|---|---|---|
| 空间标准差 | $\sigma_s$ | 2.0 | — |
| 颜色标准差 | $\sigma_r$ | 0.25 | 天气纹理保留 |
| 法线标准差 | $\sigma_n$ | 0.5 | 几何边缘保护 |
| 窗口半径 | $r$ | 3 | — |

### 9.2 Reinhard 色调映射

$$
L_{\text{mapped}} = \frac{L \cdot e \cdot \left(1 + \dfrac{L \cdot e}{w^2}\right)}{1 + L \cdot e}
$$

随后进行 sRGB 分段 OETF 校正：

$$
L_{\text{out}} = \begin{cases}
12.92\,L_{\text{mapped}}, & L_{\text{mapped}} \le 0.0031308 \\[4pt]
1.055\,L_{\text{mapped}}^{1/2.4} - 0.055, & \text{否则}
\end{cases}
$$

### 9.3 多级高斯 Bloom

提取亮部：

$$
B_0 = \max(I - t,\ 0)
$$

三级高斯模糊，标准差 $\sigma_k = 1.5 \cdot 2^{\,k-1}$，权重 $w_k = 1/2^k$（$k = 1,2,3$）：

$$
I_{\text{final}} = \text{clamp}\!\left(I + s \sum_{k=1}^{3} w_k \cdot \text{GaussianBlur}_{\sigma_k}(B_0),\ 0,\ 1\right)
$$

### 9.4 参数化场景配置（Day / Night）

| 参数 | 白天 (day) | 夜晚 (night) |
|---|---|---|
| Bloom 阈值 $t$ | 0.75 | 0.15 |
| Bloom 强度 $s$ | 0.05 | 0.35 |
| 曝光度 $e$ | 1.2 | 1.8 |
| 白点 $w$ | 5.0 | 3.0 |
| 噪声基准 $\sigma$ | 0.05 | 0.08 |
| 天气粒子 | 关闭 | 雨雪开启 |

主循环通过统一配置表查表，零硬编码分支。

---

## 10. 评估指标

全参考评价采用亮度通道计算。

**SSIM**（结构相似性），$7\times7$ 均匀滤波核：

$$
\text{SSIM} = \frac{(2\mu_x\mu_y + C_1)(2\sigma_{xy} + C_2)}{(\mu_x^2 + \mu_y^2 + C_1)(\sigma_x^2 + \sigma_y^2 + C_2)},\quad C_1=(0.01)^2,\ C_2=(0.03)^2
$$

**MSE**：

$$
\text{MSE} = \frac{1}{N}\sum (I - I_{\text{gt}})^2
$$

**梯度幅值**：

$$
\|\nabla L\| = \sqrt{\left(\frac{\partial L}{\partial x}\right)^2 + \left(\frac{\partial L}{\partial y}\right)^2}
$$

---

## 11. 附录：符号表与默认参数

### 11.1 主要符号

| 符号 | 定义 | 量纲 |
|---|---|---|
| $G$ | 顶点池 | 混合 |
| $T$ | 变换矩阵表 | 混合 |
| $I$ | 全局索引缓冲 | [1] |
| $M_j$ | 4×4 齐次矩阵 | 混合 |
| $k(N)$ | 自适应节点宽度 | [1] |
| $\sigma(N)$ | 法线最大偏离角 | [rad] |
| $V_{\mathbf{g}}$ | 局部团顶点集 | [L] |
| $\mathbf{c}_{\mathbf{g}}$ | 局部团形心 | [L] |
| $\mathbf{o}$ | 全局顶点均值 | [L] |
| $\mathbf{c}_{\text{cam}}$ | 相机位置 | [L] |
| $E$ | 曝光因子 | [1] |
| $A(\mathbf{u})$ | Albedo Buffer | [1] |

### 11.2 DyVeRT v1.4 默认参数表

| 参数 | 值 | 说明 |
|---|---|---|
| $k_{\min}, k_{\max}$ | 2, 8 | BVH 宽度范围 |
| $\gamma$ (宽度映射) | 1.5 | 陡峭度 |
| $\sigma_{\max}$ | $\pi/2$ | 最大法线偏离角 |
| 叶子容量 | 6 | 局部 BVH 叶三角形数 |
| $h$ 格宽系数 | $2\sim3$ 倍平均边长 | 相机空间分块 |
| $\epsilon$ 体积阈值 | $0.05h$ | 封闭检验 |
| $\lambda$ 外推系数 | $1.5\sim2.5$ | 相机外推 |
| $\sigma_s$ | 2.0 | 滤波空间标准差 |
| $\sigma_r$ | 0.25 | 滤波颜色标准差 |
| $\sigma_n$ | 0.5 | 滤波法线标准差 |
| $r$ | 3 | 滤波窗口半径 |
| **Gamma 校正方式** | **sRGB 分段 OETF** | 暗部线性 + 亮部幂函数 |

---

*项目名称：DyVeRT (Dynamic Vectorized Ray Tracer)*  
*文档版本：v1.4*  
*状态：全采样渲染 + 跨域双边滤波管线锁定，可交付基线*

