# 通用辛积分器框架下的 HMC 算法

## 1. 哈密顿体系
目标分布 \(\pi(x) \propto e^{-S(x)}\)，构造哈密顿量

\[
H(x, p) = S(x) + \frac{1}{2} p^{\!\top} M^{-1} p,
\]

其中 \(M\) 为正定质量矩阵（通常取 \(M = I\)）。力定义为作用量的负梯度

\[
F(x) = -\nabla S(x).
\]

正则方程

\[
\dot x = M^{-1} p, \qquad \dot p = F(x)
\]

产生保积（辛）流。

## 2. 基本辛映射
- **位置推进** \(\Phi_h^{\rm q}\)：冻结动量，沿群作用更新位置  
  \[
  \Phi_h^{\rm q} : (x, p) \mapsto \Bigl( \mathcal P\bigl[ x + h\,M^{-1} p \bigr],\; p \Bigr).
  \]
  \(\mathcal P\) 将自变量投影到 \(x\) 所在的约束流形；若无约束则为恒等。

- **动量推进** \(\Psi_h^{\rm p}\)：冻结位置，沿力方向更新动量  
  \[
  \Psi_h^{\rm p} : (x, p) \mapsto \Bigl( x,\; p + h\,F(x) \Bigr).
  \]

## 3. 蛙跳积分器
取半步动量包围整步位置：

\[
\mathcal T = \Psi_{\Delta t/2}^{\rm p}
            \circ \bigl( \Phi_{\Delta t}^{\rm q} \circ \Psi_{\Delta t}^{\rm p} \bigr)^{L-1}
            \circ \Phi_{\Delta t}^{\rm q}
            \circ \Psi_{\Delta t/2}^{\rm p}.
\]

等价地，若按部分动量拆分半步映射 \(\Theta_h\)，可写为

\[
\mathcal T = \Theta_{\Delta t/2}
            \circ \bigl( \Psi_{\Delta t}^{\rm p} \circ \Phi_{\Delta t}^{\rm q} \bigr)^{L}
            \circ \Theta_{\Delta t/2}.
\]

## 4. Metropolis 接受准则
对初态 \((x_0, p_0)\) 执行轨迹映射 \((x_1, p_1) = \mathcal T(x_0, p_0)\)，计算能量差

\[
\Delta H = H(x_1, p_1) - H(x_0, p_0),
\]

接受概率

\[
\alpha = \min\!\bigl(1,\; e^{-\Delta H}\bigr).
\]

## 5. 规范场 \((U, \gamma)\) 的具体实现

### 变量与动能
- \(U \in SU(N)\)，共轭动量 \(P_U \in \mathfrak{su}(N)\)（迹零厄米矩阵）；
- \(\gamma \in \mathbb R^d\)，共轭动量 \(P_\gamma \in \mathbb R^d\)；
- 动能  
  \[
  T = \frac12 \|P_U\|_F^2 + \frac12 \|P_\gamma\|_2^2.
  \]

### 位置推进 \(\Phi_\varepsilon^{\rm q}\)
\[
\Phi_\varepsilon^{\rm q}(U, \gamma, P_U, P_\gamma) =
\Bigl(
  \mathcal P_{SU(N)}\!\bigl[ U e^{i\varepsilon P_U} \bigr],\;
  \mathcal W[\gamma + \varepsilon P_\gamma],\;
  P_U,\; P_\gamma
\Bigr),
\]
其中
- \(e^{i\varepsilon P_U}\) 为矩阵指数（\(iP_U\) 为 \(\mathfrak{su}(N)\) 元素）；
- \(\mathcal P_{SU(N)}\) 为到 \(SU(N)\) 的最近点投影（通常由极分解 \(W = U\Sigma V^\dagger \mapsto UV^\dagger / \det^{1/N}\) 实现）；
- \(\mathcal W\) 为 \(\gamma\) 所在流形的投影，无约束时取恒等。

### 动量推进 \(\Psi_\varepsilon^{\rm p}\)
\[
\Psi_\varepsilon^{\rm p}(U, \gamma, P_U, P_\gamma) =
\Bigl(
  U,\; \gamma,\;
  P_U + \varepsilon \widetilde F_U(x),\;
  P_\gamma + \varepsilon F_\gamma(x)
\Bigr),
\]
其中
- \(F_\gamma(x) = -\nabla_\gamma S\) 为标量场力；
- \(F_U(x)\) 为规范场力，由左平移生成元定义：对任意 \(P_U \in \mathfrak{su}(N)\)，
  \[
  \frac{\mathrm d}{\mathrm d\varepsilon} S\!\bigl( U e^{i\varepsilon P_U}, \gamma \bigr) \Big|_{\varepsilon=0}
  = \operatorname{tr}\!\bigl( F_U \, P_U \bigr);
  \]
- \(\widetilde F_U = F_U + F_{\rm extra}\)，\(F_{\rm extra}\) 为可选附加力（默认零）。

### 示例：SU(2) 上的 von Mises–Fisher 型作用量
若 \(S_U(U) = -\beta\,\operatorname{Re}\operatorname{tr}(U_{\rm target}^\dagger U)\)，则

\[
F_U = -\beta\,\operatorname{Proj}_{\mathfrak{su}(2)}\!\Bigl[
        \operatorname{Im}\bigl( U_{\rm target}^\dagger U \bigr)
      \Bigr],
\]

其中 \(\operatorname{Proj}_{\mathfrak{su}(2)}[A] = A - \frac12 \operatorname{tr}(A)\,I\) 为到迹零厄米矩阵的投影。

---

组合上述映射即得轨迹映射 \(\mathcal T\)，代入第 4 节的接受准则完成 HMC 单步更新。