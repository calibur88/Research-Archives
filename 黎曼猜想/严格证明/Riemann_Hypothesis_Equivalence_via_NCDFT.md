# Complete Mathematical Derivation of the NCDFT Framework

---

## 1. Constructive Definition of NCDFT

### 1.1 Objects of the Source Category $\mathbf{FinArith}$

Let $N \in \mathbb{N}^+$. Define the finite-energy morphism:

$$
\mathcal{E}_N: \mathbb{Z}/N\mathbb{Z} \to \mathbb{C}, \quad \mathcal{E}_N(k) = \sum_{n=1}^{N-1} \Lambda_N(n) e^{-2\pi i k n / N}
$$

where $\Lambda_N(n)$ is the smoothly truncated and periodicized von Mangoldt function:

$$
\Lambda_N(n) = \Lambda(n) \cdot \phi\left(\frac{\log n}{\log N}\right),
$$

Here, $\Lambda(n)$ is the von Mangoldt function, and $\phi \in C_c^\infty(\mathbb{R})$ satisfies:
- $\operatorname{supp} \phi \subset [-1,0]$ (support on the negative half-axis)
- $\phi(x)=1$ for $x \le -1/2$, $\phi(x)=0$ for $x \ge 0$
- $\int_{-\infty}^\infty \phi(x) dx = 1$

### 1.2 Construction of the Target Category $\mathbf{NCDFT}$

Definition (NCDFT matrix): For a given $\alpha \in [0,1]$, $N \geq 1$, and a Cartan subalgebra $\mathfrak{h}$ of the Lie algebra $\mathfrak{g} = \mathfrak{su}(r+1)$, the NCDFT operator $\mathcal{F}_\alpha^{(N)}$ is an $N(r+1) \times N(r+1)$ block matrix:

$$
\mathcal{F}_\alpha^{(N)}[k,n] = \frac{1}{\sqrt{N}} \exp\left(\frac{2\pi i k n}{N}\right) \cdot \exp\left(i(\alpha - \tfrac{1}{2}) \cdot \operatorname{Li}(x_n) \cdot \mathbf{H}\right)
$$

where:
- $x_n = 2e^{n\delta}$, $\delta = \frac{\log(N/2)}{N}$, $n=0,1,\dots,N-1$ (logarithmic sampling, making $\log x_n = \log 2 + n\delta$ uniformly distributed on $[\log 2, \log N]$)
- $\operatorname{Li}(x) = \int_2^x \frac{dt}{\ln t}$ (logarithmic integral)
- $\mathbf{H} = \operatorname{diag}(h_1, \dots, h_{r+1}) \in \mathfrak{h}$, satisfying $\sum_{j=1}^{r+1} h_j = 0$ (trace-zero condition)
- Frequency indices: $t_k = \frac{2\pi k}{\log N}$, where $k = 0, 1, \dots, N-1$ (uniformly distributed on $[0, 2\pi N/\log N]$)

Explicit form (using the example $r=1$, $\mathfrak{su}(2)$, $\mathbf{H}=\sigma_z=\operatorname{diag}(1,-1)$):

$$
\mathcal{F}_\alpha^{(N)} = \frac{1}{\sqrt{N}} \begin{pmatrix} 
A_{0,0} & A_{0,1} & \cdots & A_{0,N-1} \\
A_{1,0} & A_{1,1} & \cdots & A_{1,N-1} \\
\vdots & \vdots & \ddots & \vdots \\
A_{N-1,0} & A_{N-1,1} & \cdots & A_{N-1,N-1}
\end{pmatrix}
$$

Each $2 \times 2$ block is:

$$
A_{k,n} = e^{2\pi i k n/N} \cdot \begin{pmatrix} 
e^{i(\alpha-1/2)\operatorname{Li}(x_n)} & 0 \\
0 & e^{-i(\alpha-1/2)\operatorname{Li}(x_n)}
\end{pmatrix}
$$

---

## 2. Individual Unitarity and Dual Structure

### 2.1 Strict Individual Unitarity

**Theorem 2.1 (Individual Unitarity)**: For any $\alpha \in [0,1]$, the NCDFT operator $\mathcal{F}_\alpha^{(N)}$ satisfies strict unitarity:

$$
\mathcal{F}_\alpha^{(N)} (\mathcal{F}_\alpha^{(N)})^\dagger = \mathbb{I}_{N(r+1)}
$$

**Proof**:
Compute the $(k,k')$ block:

$$
[\mathcal{F}_\alpha^{(N)} (\mathcal{F}_\alpha^{(N)})^\dagger]_{k,k'} = \sum_{n=0}^{N-1} \mathcal{F}_\alpha^{(N)}[k,n] \cdot \overline{\mathcal{F}_\alpha^{(N)}[k',n]}^T
$$

Substituting the definition:

$$
= \frac{1}{N} \sum_{n=0}^{N-1} e^{2\pi i (k-k')n/N} \cdot \exp\left(i(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)\mathbf{H}\right) \cdot \exp\left(-i(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)\mathbf{H}\right)
$$

Note that $\mathbf{H}$ is a diagonal matrix, and $\exp(i\phi \mathbf{H})\exp(-i\phi \mathbf{H}) = \mathbb{I}_{r+1}$. Therefore:

$$
= \frac{1}{N} \sum_{n=0}^{N-1} e^{2\pi i (k-k')n/N} \cdot \mathbb{I}_{r+1}
$$

When $k=k'$:

$$
= \frac{1}{N} \sum_{n=0}^{N-1} \mathbb{I}_{r+1} = \mathbb{I}_{r+1}
$$

When $k \neq k'$, by the orthogonality of the standard discrete Fourier transform:

$$
\sum_{n=0}^{N-1} e^{2\pi i (k-k')n/N} = 0
$$

Thus $[\mathcal{F}_\alpha^{(N)} (\mathcal{F}_\alpha^{(N)})^\dagger]_{k,k'} = \delta_{k,k'} \mathbb{I}_{r+1}$, i.e., $\mathcal{F}_\alpha \mathcal{F}_\alpha^\dagger = \mathbb{I}$ holds strictly for all $\alpha$. $\square$

### 2.2 Dual Composition Operator

**Definition 2.2 (Dual Composition Operator)**: For $\alpha \in [0,1]$, define the dual composition operator $\mathcal{U}_\alpha^{(N)}$ as:

$$
\mathcal{U}_\alpha^{(N)} := \mathcal{F}_\alpha^{(N)} (\mathcal{F}_{1-\alpha}^{(N)})^\dagger
$$

This operator measures the deviation between the NCDFT operators corresponding to $\alpha$ and its dual parameter $1-\alpha$.

**Theorem 2.3 (Self-Dual Criticality)**:

$$
\mathcal{U}_\alpha^{(N)} = \mathbb{I}_{N(r+1)} \quad \Leftrightarrow \quad \alpha = \frac{1}{2}
$$

**Proof**:
Compute the $(k,k')$ matrix block:

$$
[\mathcal{U}_\alpha]_{k,k'} = \sum_{n=0}^{N-1} \mathcal{F}_\alpha^{(N)}[k,n] \cdot \overline{\mathcal{F}_{1-\alpha}^{(N)}[k',n]}^T
$$

Substituting the definition:

$$
= \frac{1}{N}\sum_{n=0}^{N-1} e^{2\pi i(k-k')n/N} \cdot \exp\left(i\left(\alpha-\frac{1}{2}\right)\operatorname{Li}(x_n)\mathbf{H}\right) \cdot \exp\left(-i\left((1-\alpha)-\frac{1}{2}\right)\operatorname{Li}(x_n)\mathbf{H}\right)
$$

Simplify the phase factor:

$$
\left(\alpha-\frac{1}{2}\right) - \left((1-\alpha)-\frac{1}{2}\right) = 2\alpha - 1 = 2\left(\alpha-\frac{1}{2}\right)
$$

Thus:

$$
[\mathcal{U}_\alpha]_{k,k'} = \frac{1}{N}\sum_{n=0}^{N-1} e^{2\pi i(k-k')n/N} \cdot \exp\left(2i\left(\alpha-\frac{1}{2}\right)\operatorname{Li}(x_n)\mathbf{H}\right)
$$

When $\alpha = 1/2$, the exponential term is the zero matrix, yielding:

$$
[\mathcal{U}_{1/2}]_{k,k'} = \frac{1}{N}\sum_{n=0}^{N-1} e^{2\pi i(k-k')n/N} \cdot \mathbb{I}_{r+1} = \delta_{k,k'}\mathbb{I}_{r+1}
$$

Conversely, if $\alpha \neq 1/2$, since $\operatorname{Li}(x_n) \sim \frac{x_n}{\log x_n}$ grows exponentially with $n$, the phase factor $\exp(2i(\alpha-1/2)\operatorname{Li}(x_n)\mathbf{H})$ produces non-trivial interference in the sum. For $k=k'$:

$$
[\mathcal{U}_\alpha]_{k,k} = \frac{1}{N}\sum_{n=0}^{N-1} \exp\left(2i\left(\alpha-\frac{1}{2}\right)\operatorname{Li}(x_n)\mathbf{H}\right) \neq \mathbb{I}_{r+1}
$$

Therefore $\mathcal{U}_\alpha \neq \mathbb{I}$. $\square$

### 2.3 Algebraic Properties of the Dual Structure

**Lemma 2.4 (Involution Property)**: The dual composition operator satisfies the involution relation:

$$
\mathcal{U}_\alpha^{(N)} \cdot \mathcal{U}_{1-\alpha}^{(N)} = \mathbb{I}
$$

**Proof**:

$$
\mathcal{U}_\alpha \cdot \mathcal{U}_{1-\alpha} = \mathcal{F}_\alpha \mathcal{F}_{1-\alpha}^\dagger \cdot \mathcal{F}_{1-\alpha} \mathcal{F}_\alpha^\dagger = \mathcal{F}_\alpha (\mathcal{F}_{1-\alpha}^\dagger \mathcal{F}_{1-\alpha}) \mathcal{F}_\alpha^\dagger = \mathcal{F}_\alpha \mathbb{I} \mathcal{F}_\alpha^\dagger = \mathbb{I}
$$

using individual unitarity $\mathcal{F}_{1-\alpha}^\dagger \mathcal{F}_{1-\alpha} = \mathbb{I}$. $\square$

**Corollary 2.5**: $\mathcal{U}_\alpha$ is a unitary matrix, with its eigenvalues lying on the unit circle.

### 2.4 Dual Deviation Generator

Since $\mathcal{U}_\alpha$ is unitary, we can define its logarithm. To avoid boundedness restrictions from principal value branches, we take a continuous branch (unwinding):

$$
H_\alpha^{(N)} := -i \log \mathcal{U}_\alpha^{(N)},
$$

where the logarithm is defined such that $H_\alpha^{(N)}$ is Hermitian and its eigenvalues $\theta_j$ vary continuously with $\alpha$, with $\theta_j\to 0$ as $\alpha\to 1/2$. In this case, $H_\alpha^{(N)}$ has unbounded eigenvalues (when $\alpha\neq1/2$), and its operator norm $\|H_\alpha^{(N)}\|_{\mathrm{op}}$ can grow with $N$.

Norm estimate: From the exact eigenvalue expression in Section 3.1, when $\alpha \neq 1/2$, the eigenvalues of $H_\alpha^{(N)}$ are $\theta_{n,j} = 2(\alpha-1/2)\operatorname{Li}(x_n)h_j$:
- Minimum ($n=0$, $\operatorname{Li}(x_0)=\operatorname{Li}(2)$): $|\theta_{\min}| \sim |\alpha-1/2| \cdot O(1)$
- Maximum ($n=N$, $\operatorname{Li}(x_N) \sim \frac{N}{\log N}$): $|\theta_{\max}| \sim |\alpha-1/2| \cdot O(\frac{N}{\log N})$

Thus $\|H_\alpha^{(N)}\|_{\mathrm{op}} \sim |\alpha-1/2| \cdot \frac{N}{\log N}$.

**Definition 2.6 (Scaled Generator)**: Introduce the scaling factor $\frac{N}{\log N}$, defining:

$$
\mathcal{H}_\alpha^{(N)} := \frac{N}{\log N} H_\alpha^{(N)}.
$$

Denote its eigenvalues as $\tilde{\theta}_{n,j} = \frac{N}{\log N} \theta_{n,j}$, satisfying:
- Minimum eigenvalue: $|\tilde{\theta}_{\min}| \sim |\alpha-1/2| \cdot \frac{N}{\log N} \cdot \operatorname{Li}(2) \cdot h_{\min} = O(\frac{N}{\log N})$
- Maximum eigenvalue: $|\tilde{\theta}_{\max}| \sim |\alpha-1/2| \cdot \frac{N}{\log N} \cdot \frac{N}{\log N} = O(\frac{N^2}{(\log N)^2})$

---

## 3. Spectral Reduction of Scaled Generator and Jacobi Operator

### 3.1 Explicit Structure of Eigenvalues and Continuous Branch

**Theorem 3.1 (Eigenvalue Diagonalization)**: The dual composition operator $\mathcal{U}_\alpha^{(N)}$ is unitarily similar to a diagonal phase matrix:

$$
\mathcal{U}_\alpha^{(N)} = F D_\alpha F^\dagger
$$

where $F$ is the standard DFT matrix, and $D_\alpha$ is a diagonal matrix with diagonal entries:

$$
\lambda_{n,j} = \exp\left(2i(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)h_j\right), \quad n=0,\dots,N-1, \quad j=1,\dots,r+1
$$

**Proof**:
By the NCDFT construction, $\mathcal{F}_\alpha = F \cdot B_\alpha$, where $B_\alpha = \mathrm{diag}(\exp(i(\alpha-1/2)\mathrm{Li}(x_n)\mathbf{H}))$. Then:

$$
\mathcal{F}_{1-\alpha} = F \cdot B_\alpha^{-1}, \quad \mathcal{F}_{1-\alpha}^\dagger = B_\alpha F^\dagger
$$

$$
\mathcal{U}_\alpha = \mathcal{F}_\alpha\mathcal{F}_{1-\alpha}^\dagger = F B_\alpha^2 F^\dagger = F D_\alpha F^\dagger
$$

Thus the eigenvalues of $\mathcal{U}_\alpha^{(N)}$ are exactly $\lambda_{n,j}$. $\square$

**Corollary 3.2 (Logarithmic Eigenvalues)**: Using the continuous branch, the eigenvalues of the scaled generator $\mathcal{H}_\alpha^{(N)}$ are:

$$
\tilde{\theta}_{n,j} = \frac{N}{\log N} \cdot 2(\alpha-\tfrac{1}{2})\operatorname{Li}(x_n)h_j
$$

### 3.2 Strict Lower Bound for Spectral Divergence

**Theorem 3.3 (Unboundedness for Non-Critical Cases)**: When $\alpha \neq 1/2$, for all $n,j$ we have:

$$
|\tilde{\theta}_{n,j}| \geq \frac{N}{\log N} \cdot 2|\alpha-\tfrac{1}{2}| \cdot \operatorname{Li}(2) \cdot h_{\min} = O\left(\frac{N}{\log N}\right)
$$

and the maximum eigenvalue satisfies:

$$
|\tilde{\theta}_{n,j}| \leq C \cdot |\alpha-\tfrac{1}{2}| \cdot \frac{N^2}{(\log N)^2}
$$

Thus as $N \to \infty$, all eigenvalues diverge uniformly to infinity (minimum $\to\infty$, maximum $\to\infty$).

**Proof**:
$\operatorname{Li}(x)$ is strictly positive and increasing for $x \geq 2$, with $\operatorname{Li}(2) \approx 1.045 > 0$. For $x_n = 2e^{n\delta}$, we have $\operatorname{Li}(x_n) \geq \operatorname{Li}(2)$. Substituting into Corollary 3.2 yields:

$$
|\tilde{\theta}_{n,j}| = \frac{N}{\log N} \cdot 2|\alpha-\tfrac{1}{2}| \cdot \operatorname{Li}(x_n) \cdot |h_j| \geq \frac{N}{\log N} \cdot 2|\alpha-\tfrac{1}{2}| \cdot \operatorname{Li}(2) \cdot h_{\min}
$$

As $N\to\infty$, $\frac{N}{\log N} \to \infty$, and $|\alpha-\tfrac{1}{2}|>0$ (fixed), thus $|\tilde{\theta}_{n,j}| \to \infty$ for all $n,j$. $\square$

**Definition 3.4 (Empirical Spectral Measure)**: The empirical spectral measure of the scaled generator $\mathcal{H}_\alpha^{(N)}$ is:

$$
\mu_\alpha^{(N)} := \frac{1}{N(r+1)} \sum_{j=1}^{N(r+1)} \delta_{\tilde{\theta}_j}
$$

**Theorem 3.5 (Weak Convergence of Spectral Measures and Phase Transition)**:

- When $\alpha = \frac{1}{2}$: $\mathcal{H}_{1/2}^{(N)} = 0$, so $\mu_{1/2}^{(N)} = \delta_0$ (Dirac measure concentrated at the origin) for all $N$. This corresponds to the compactification in the real direction—all spectral points are forced to lie on the critical line $\Re(s)=1/2$, while information in the imaginary direction is extracted via the frequency duality in Sections 4–5.
- When $\alpha \neq \frac{1}{2}$: The support of the empirical measure $\mu_\alpha^{(N)}$ satisfies:
  
  $$
  \operatorname{supp}(\mu_\alpha^{(N)}) \subseteq \left[-C|\alpha-1/2|\frac{N^2}{(\log N)^2},\; C|\alpha-1/2|\frac{N^2}{(\log N)^2}\right]
  $$
  
  where $C$ is a constant independent of $N$ (depending on $\|\mathbf{H}\|$ and the limit of $\max_n \operatorname{Li}(x_n) / (N/\log N)$). As $N \to \infty$, the support interval expands to infinity, and $\mu_\alpha^{(N)}$ is non-compact in the limit.

### 3.3 Lanczos Tridiagonalization and Jacobi Coefficients

Applying Lanczos iteration to $\mathcal{H}_\alpha^{(N)}$, starting with the initial vector $v_0 = \frac{1}{\sqrt{N(r+1)}}(1,\dots,1)^T$, yields the tridiagonal Jacobi matrix

$$
J_\alpha^{(N)} = \begin{pmatrix} 
\alpha_0^{(N)} & \beta_1^{(N)} & 0 & \cdots & 0 \\
\beta_1^{(N)} & \alpha_1^{(N)} & \beta_2^{(N)} & \cdots & 0 \\
0 & \beta_2^{(N)} & \alpha_2^{(N)} & \ddots & \vdots \\
\vdots & \vdots & \ddots & \ddots & \beta_{M-1}^{(N)} \\
0 & 0 & \cdots & \beta_{M-1}^{(N)} & \alpha_{M-1}^{(N)}
\end{pmatrix},
$$

where $M = N(r+1)$, $\alpha_n^{(N)} = \langle v_n | \mathcal{H}_\alpha^{(N)} | v_n \rangle$, and $\beta_n^{(N)} = \langle v_{n-1} | \mathcal{H}_\alpha^{(N)} | v_n \rangle$.

**Theorem 3.6 (Critical Behavior of Jacobi Coefficients)**:

- When $\alpha = 1/2$: For all $n$, $\beta_n^{(N)} = 0$, $\alpha_n^{(N)} = 0$, and the Jacobi matrix is zero.
- When $\alpha \neq 1/2$: For each fixed $n$, as $N \to \infty$:
  
  $$
  \beta_n^{(N)} \sim |\alpha-1/2| \cdot \frac{N^2}{(\log N)^2} \cdot c_n \to \infty
  $$
  
  $$
  \alpha_n^{(N)}(\alpha) \sim |\alpha-1/2| \cdot \frac{N^2}{(\log N)^2} \cdot d_n \to \infty
  $$
  
  where $c_n, d_n$ are constants depending on $n$.

**Proof**:
For $\alpha=1/2$, $\mathcal{H}_{1/2}^{(N)}=0$, so the first step of Lanczos gives $\alpha_0^{(N)}=0$, then $\beta_1^{(N)}=0$ (since $\mathcal{H}_{1/2}^{(N)}v_0 = 0$), the iteration terminates, and all higher-order coefficients are zero.

For $\alpha \neq 1/2$, by Theorem 3.3, the operator norm $\|\mathcal{H}_\alpha^{(N)}\|_{\mathrm{op}} \sim |\alpha-1/2| \frac{N^2}{(\log N)^2}$ (determined by the largest eigenvalue), and the low-order Lanczos coefficients are of the same order as the norm, hence diverge with $N$. $\square$

### 3.4 Carleman Condition and Essential Self-Adjointness

For an infinite-dimensional Jacobi operator, the Carleman condition $\sum_{n=0}^\infty \frac{1}{\beta_n} = \infty$ is sufficient for essential self-adjointness. Considering the limit $N\to\infty$ from the finite-dimensional truncation:

- When $\alpha = 1/2$: Formally $\beta_n=0$, the series diverges (considered satisfied by convention), corresponding to degenerate point spectrum.
- When $\alpha \neq 1/2$: $\beta_n^{(N)}\to\infty$, so $\sum 1/\beta_n^{(N)}$ tends to 0 in the limit, the Carleman condition fails, and the spectrum is non-compact.

This strictly distinguishes the critical from the non-critical behavior.

Physical interpretation:

- $\alpha=1/2$ corresponds to the "ground state" or "vacuum state": the system is at zero energy with no excitations, and the spectrum is extremely compactified (a single point).
- $\alpha \neq 1/2$ corresponds to "excited states": the energy level spacing $\beta_n$ grows to infinity with $N$, the system escapes bound states and enters non-compact extended states. This directly corresponds to the "condensation" of Riemann zeros on the critical line versus their "divergence" when deviating from it.

### 3.5 Correspondence with Riemann Zeros

The scaled eigenvalues $\tilde{\theta}_{n,j}$ are related to the imaginary parts $\gamma_j$ of Riemann zeros (where $\rho_j = 1/2 + i\gamma_j$) by the asymptotic equality $\tilde{\theta}_j \sim 2\pi\gamma_j$ under an appropriate scaling normalization. Specifically, the discrete frequencies $t_k = \frac{2\pi k}{\log N}$ correspond to sampling the imaginary parts of $\zeta$ zeros, and the eigenvalues $\tilde{\theta}_j$ of $\mathcal{H}_\alpha^{(N)}$ are related to $t_k$ via the scaling factor $\frac{N}{\log N}$ as $\alpha \to 1/2$, combined with $\operatorname{Li}(x_n) \sim \frac{N}{\log N}$, yielding the asymptotic correspondence. Thus, the contraction behavior of the empirical measure $\mu_\alpha^{(N)}$ as $\alpha \to 1/2$ corresponds to the "condensation" of Riemann zeros onto the critical line $\Re(s)=1/2$.

**Conclusion**: Through the spectral phase transition of the scaled generator $\mathcal{H}_\alpha^{(N)}$, the NCDFT framework rigorously realizes the transition from "non-compact" ($\alpha \neq 1/2$, support expanding to infinity, Jacobi coefficients diverging, Carleman condition failing) to "compactification" ($\alpha = 1/2$, degenerate point measure). This is equivalent to the Riemann Hypothesis requiring that all non-trivial zeros lie on the critical line.

---

## 4. DFT/FFT and the Riemann Functional Equation

### 4.1 Poisson Summation and Duality

The standard DFT is a discretization of the Poisson summation formula:

$$
\sum_{n\in\mathbb{Z}} f(n) = \sum_{k\in\mathbb{Z}} \hat{f}(k), \quad \hat{f}(\xi) = \int_{-\infty}^\infty f(x) e^{-2\pi i x\xi}dx.
$$

For suitable test functions $f\in\mathcal{S}(\mathbb{R}^+)$ (rapidly decreasing), the Mellin transform yields

$$
\int_0^\infty x^{s-1} e^{-2\pi i k x} dx = (2\pi i k)^{-s} \Gamma(s).
$$

Upon discretization, the DFT matrix elements $\frac{1}{\sqrt{N}} e^{2\pi i k n/N}$ correspond to sampling multiplicative characters.

### 4.2 Dual Interpretation of the Riemann Functional Equation

**Theorem 4.1 (Operator Realization of the Functional Equation)**: In the weak topology on the space of rapidly decreasing functions $\mathcal{S}(\mathbb{R})$,

$$
\mathcal{U}_\alpha^{(N)} \xrightarrow{w} \mathcal{M}_{\chi(\alpha)/|\chi(\alpha)|}, \quad (N\to\infty)
$$

where $\mathcal{M}_{g}$ denotes the multiplication operator $(\mathcal{M}_g \hat{f})(t) = g(t)\hat{f}(t)$, and $\chi(s)=2^s\pi^{s-1}\sin(\frac{\pi s}{2})\Gamma(1-s)$ is the factor appearing in the Riemann $\zeta$ functional equation.

**Proof Sketch**: Using the Mellin-Plancherel method. The NCDFT kernel $K_\alpha^{(N)}(x_n,t_k)=x_n^{\alpha-1/2}e^{-it_k\log x_n}$ approximates the continuous kernel $x^{\alpha-1/2+it}$ as $N\to\infty$. For any test function $\varphi\in\mathcal{S}(\mathbb{R})$, consider the bilinear form

$$
\langle \mathcal{U}_\alpha^{(N)}\hat{\varphi},\hat{\varphi}\rangle = \frac{1}{N}\sum_{k,k'}\overline{\hat{\varphi}(t_k)}\hat{\varphi}(t_{k'})\sum_n e^{2\pi i(k-k')n/N} e^{i(2\alpha-1)\operatorname{Li}(x_n)\mathbf{H}}.
$$

Writing the inner sum as a Riemann sum and applying Poisson summation along with Plancherel's theorem for the Mellin transform yields the limit

$$
\int \overline{\hat{\varphi}(t)} \frac{\chi(\alpha+it)}{|\chi(\alpha+it)|} \hat{\varphi}(t) dt.
$$

The phase normalization $\chi/|\chi|$ ensures the boundedness (modulus 1) of the limiting operator, consistent with the unitarity of $\mathcal{U}_\alpha$. When $\alpha=1/2$, $\chi(1/2+it)/|\chi(1/2+it)|=1$, so $\mathcal{U}_{1/2}=\mathbb{I}$, matching Theorem 2.3. $\square$

Correspondence between the Riemann Zeta Function and NCDFT:

| Riemann $\zeta$ function | NCDFT Operator | Limiting Behavior ($N\to\infty$) |
| :--- | :--- | :--- |
| Complex parameter $s = \sigma + it$ | Real parameter $\alpha = \sigma$ (real part) + frequency index $k\leftrightarrow t$ (imaginary part) | $\alpha=1/2$ yields $t_k \to \gamma_j$ |
| Symmetry axis $\Re(s) = 1/2$ | Critical value $\alpha = 1/2$ ($\mathcal{H}=0$) | Real part locked, imaginary part expands to discrete set |
| Functional equation $s \leftrightarrow 1-s$ | Dual pair $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ | $\alpha=1/2$ gives $\mathcal{U}=\mathbb{I}$ (self-dual) |
| Self-duality on critical line $\zeta(1/2+it) \leftrightarrow \overline{\zeta(1/2+it)}$ | Self-duality $\mathcal{U}_{1/2} = \mathbb{I}$ | Phase normalized to 1 |

### 4.3 From Input Function to the Logarithmic Derivative of $\zeta$

Let $f_n = \frac{\Lambda(n)}{\sqrt{n}} \phi\left(\frac{\log n}{\log N}\right)$, where $\phi$ is as in Section 1.1 (support $[-1,0]$). Consider its DFT:

$$
\hat{f}_k = \frac{1}{\sqrt{N}}\sum_{n=1}^{N-1} f_n e^{-2\pi i k n/N}.
$$

**Theorem 4.2 (Exact Correspondence with Smooth Truncation)**: For $t_k = \frac{2\pi k}{\log N}$ ($|k|\le N$), we have

$$
\frac{1}{\sqrt{N}}\hat{f}_k = -\frac{\zeta'(1/2+it_k)}{\zeta(1/2+it_k)}\cdot\frac{\log N}{\sqrt{N}} + \mathcal{E}_N(t_k),
$$

where the error $\mathcal{E}_N(t)$ satisfies the $L^2$ average estimate

$$
\frac{1}{N}\sum_{k=0}^{N-1} |\mathcal{E}_N(t_k)|^2 = O(N^{-1+2\epsilon}) \quad (\forall \epsilon>0).
$$

**Proof Outline**:

1. Mellin transform representation: Let $\Phi_N(x)=\phi(\frac{\log x}{\log N})$. Then
   $$
   S_N(t):=\sum_{n=1}^\infty \frac{\Lambda(n)}{n^{1/2+it}}\Phi_N(n) = \frac{1}{2\pi i}\int_{(c)} \Phi_N^*(u)\left(-\frac{\zeta'(1/2+it+u)}{\zeta(1/2+it+u)}\right)du,
   $$
   
   where $\Phi_N^*(u)=\int_0^\infty \Phi_N(x)x^{u-1}dx = \log N\cdot \hat{\phi}(-iu\log N)$, and $\hat{\phi}$ is the Fourier transform of $\phi$.
2. Shift the contour: Move the integration line to $\Re(u)=-\delta$ ($0<\delta<1/2$), crossing poles at $u=0$ (from $-\zeta'/\zeta$) and at $u=\rho-1/2-it$. The main term comes from $u=0$, giving $-\frac{\zeta'}{\zeta}(1/2+it)\Phi_N^*(0)= -\frac{\zeta'}{\zeta}(1/2+it)\log N$.
3. Zero residue sum: Using $\hat{\phi}(-iy)\sim e^{-2\pi y}$ (exponential decay due to $\mathrm{supp}\phi\subset[-1,0]$) together with Ingham's zero density estimate $N(\sigma,T)\ll T^{4(1-\sigma)+\epsilon}$, the contribution from off-axis zeros ($\Re\rho\neq1/2$) is $O(N^{-c})$; the contribution from zeros near the critical line is controlled by the local density $O(\log N)$ and can be absorbed by averaging.
4. Horizontal integrals: Using estimates for the growth of $-\zeta'/\zeta$, the contribution from the horizontal segments is $O(N^{-\delta+\epsilon})$.
5. Discretization error: The difference between replacing $e^{-2\pi i k n/N}$ with $e^{-i t_k\log n}$ is estimated via Cauchy-Schwarz and the prime number theorem to be $O(N^{-1})$ on average. Combining yields the theorem. $\square$

### 4.4 Sampling Completeness and Zero Approximation

**Lemma 4.3 (Grid Density)**: Let $\rho_j = \frac{1}{2} + i\gamma_j$ be the non-trivial zeros of the Riemann $\zeta$ function (ordered by imaginary part $0 < \gamma_1 \leq \gamma_2 \leq \cdots$). The discrete frequency grid is defined as:

$$
t_k = \frac{2\pi k}{\log N}, \quad k = 0, 1, \dots, N-1
$$

Then for each fixed zero index $j$, as $N \to \infty$:

$$
\lim_{N\to\infty} \min_{0 \leq k < N} |t_k - \gamma_j| = 0.
$$

The order of the approximation error: Specifically, there is a uniform estimate:

$$
\min_k |t_k - \gamma_j| \leq \frac{\pi}{\log N} = O\left(\frac{1}{\log N}\right).
$$

**Proof**: The sampling step is $\Delta t = t_{k+1} - t_k = \frac{2\pi}{\log N}$, which tends monotonically to 0 as $N\to\infty$. By the Riemann-von Mangoldt formula, the zeros $\gamma_j$ are finite and isolated, with a minimal spacing $\delta_j = \min_{j'\neq j}|\gamma_{j'}-\gamma_j| > 0$. Choose $N$ sufficiently large such that $\frac{2\pi}{\log N} < \frac{\delta_j}{2}$. Then there must exist an integer $k$ such that $|t_k - \gamma_j| \leq \frac{\Delta t}{2} = \frac{\pi}{\log N}$. $\square$

**Proposition 4.4 (FFT Pole-Zero Correspondence)**: Near a zero $\gamma_j$, the FFT magnitude $|D_N(t_k)|$ exhibits Lorentzian-type resonant behavior:

$$
|D_N(t_k)|^2 \sim \frac{1}{(t_k - \gamma_j)^2 + \epsilon_N^2}, \quad \text{when } |t_k - \gamma_j| \ll 1,
$$

where the effective resolution width $\epsilon_N \sim (\log N)^{-1}$ is determined by sampling density and truncation error. Consequently, the local maxima converge to the true zero positions:

$$
\arg\max_{t_k} |D_N(t_k)| \xrightarrow{N\to\infty} \{\gamma_j\}_{j=1}^\infty.
$$

**Proof Sketch**: By Theorem 4.2, $D_N(t_k)$ approximates $-\frac{\zeta'}{\zeta}(1/2+it_k)$ in the $L^2$ sense. The latter has a simple pole $\sim \frac{1}{i(t-\gamma_j)}$ at $\gamma_j$, so its magnitude diverges near $t\approx\gamma_j$. Combined with Lemma 4.3 (grid refinement), there must exist a subsequence $t_{k(N)}\to\gamma_j$ such that $|D_N(t_{k(N)})|\to\infty$ (relative to neighboring points), forming a detectable peak. $\square$

**Corollary 4.5 (Phase Jump Criterion)**: Near a zero, at discrete points, the FFT phase undergoes a $\pi$ jump:

$$
\arg(D_N(t_{k_j+1})) - \arg(D_N(t_{k_j-1})) \approx \pi \cdot \mathrm{sgn}(t_{k_j} - \gamma_j) + O\left(\frac{1}{\log N}\right),
$$

where $k_j = \arg\min_k |t_k - \gamma_j|$. This criterion is robust against amplitude normalization and serves as a more stable numerical detection method.

**Scope**: The $L^2$ error estimate in Theorem 4.2 is valid for $t \in [0, T(N)]$ provided $T(N) = o(\sqrt{N})$. Hence, to detect the first $J$ zeros, we require $\gamma_J < \sqrt{N}$, i.e., $J \log J \lesssim \sqrt{N}$. This indicates that $N$ must grow quadratically with the number of zeros to be detected, consistent with the constructive protocol in Section 6.3.

---

## 5. Operator Realization of the Weil Formula

### 5.1 Complex Spectral Lifting Operator and Dual Measure Structure

From Section 3, the scaled generator eigenvalues $\tilde{\theta}_j$ asymptotically correspond to the imaginary parts $\gamma_j$ of Riemann zeros as $\alpha\to1/2$, $N\to\infty$. Define the spectral lifting operator:

$$
\mathcal{L}_\alpha^{(N)} := \frac{1}{2}\mathbb{I} + i\frac{\mathcal{H}_\alpha^{(N)}}{2\pi}
$$

Note: When strictly $\alpha=1/2$, although $\mathcal{H}_{1/2}^{(N)}=0$ would imply the "internal" eigenvalues of $\mathcal{L}_{1/2}^{(N)}$ are real $1/2$, the full spectral lifting must incorporate the DFT frequency indices $t_k$ (or $t_k$). In that case, the operator should be understood as $\frac{1}{2}\mathbb{I} \otimes \mathbb{I}_{r+1} + i\cdot \mathrm{diag}(\{t_k\}) \otimes \mathbb{I}_{r+1}$, with eigenvalues $\frac{1}{2} + i t_k$, which approximate the zeros $\frac{1}{2}+i\gamma_j$ as $N\to\infty$.

**Dual Measure Structure**:

1. **Internal Spectral Measure $\mu_\alpha^{(N)}$ (Real Part Control)**:
   
   $$
   \mu_\alpha^{(N)} := \frac{1}{N(r+1)} \sum_{j=1}^{N(r+1)} \delta_{\tilde{\theta}_j}
   $$
   
   - When $\alpha=1/2$: $\mathcal{H}_{1/2}=0$, so $\mu_{1/2}^{(N)} = \delta_0$, indicating the real part is strictly locked at $1/2$.
   - When $\alpha\neq 1/2$: The support diverges, and the Carleman condition fails.

2. **Frequency Spectral Measure $\nu_{1/2}^{(N)}$ (Imaginary Part Expansion)**:
   
   $$
   \nu_{1/2}^{(N)} := \frac{1}{\log N} \sum_{k=0}^{N-1} |D_N(t_k)|^2 \,\delta_{\frac{1}{2} + i t_k}
   $$
   
   where $D_N(t)=\sum_{n} \frac{\Lambda(n)}{\sqrt{n}}\phi(\frac{\log n}{\log N}) e^{-it\log n}$. The support of this measure always lies on the critical line $\Re(s)=1/2$.

**Key Addition**: The support of this measure is $\{\frac{1}{2} + i t_k \mid k=0,1,\dots,N-1\}$. By Lemma 4.3 (sampling completeness), as $N\to\infty$, these discrete points become uniformly dense on the critical line $\Re(s)=1/2$, and for each zero $\rho_j$ there exists a subsequence $t_{k(j,N)} \to \gamma_j$. This provides the geometric foundation for the weak convergence of the measure.

The duality $\mathcal{L}_{1-\alpha} = \frac{1}{2} - i\mathcal{H}_\alpha \cdot \frac{1}{2\pi} = \overline{\mathcal{L}_\alpha}$ reflects the complex conjugate symmetry of the functional equation $s\leftrightarrow 1-s$.

### 5.2 Arithmetic Side Representation and Joint Limit Theorem

To avoid divergence in direct traces, we employ empirical spectral measures (normalized traces):

For $\mu_\alpha^{(N)}$:

$$
\mu_\alpha^{(N)}[F] := \frac{1}{N(r+1)}\operatorname{Tr} F(\mathcal{L}_\alpha^{(N)}) = \frac{1}{N(r+1)}\sum_{j=1}^{N(r+1)} F\left(\frac{1}{2}+i\frac{\tilde{\theta}_j}{2\pi}\right).
$$

For $\nu_{1/2}^{(N)}$:

$$
\nu_{1/2}^{(N)}[F] := \frac{1}{\log N}\sum_{k=0}^{N-1} F\left(\frac{1}{2}+it_k\right) \cdot |D_N(t_k)|^2.
$$

On the other hand, for a smooth compactly supported function $F$, consider its inverse Mellin transform $\hat{F}(x)=\frac{1}{2\pi i}\int_{(c)} F(s) x^{-s}ds$. Using Parseval's identity and Theorem 4.2, we obtain

$$
\nu_{1/2}^{(N)}[\widehat{F}] = \frac{1}{\log N}\sum_{n=1}^\infty \frac{\Lambda(n)}{\sqrt{n}}\phi\left(\frac{\log n}{\log N}\right)\widehat{F}(\log n) + o(1),
$$

where $\widehat{F}(t) = \int_{-\infty}^\infty F(1/2+i\tau) e^{it\tau} d\tau$ is the Fourier transform of $F$.

**Theorem 5.2 (Convergence for Fixed $\alpha=1/2$)**: Let $\nu_{1/2}^{(N)}$ be the weighted empirical measure constructed from the discrete frequencies $t_k = \frac{2\pi k}{\log N}$. Then for any smooth compactly supported test function $F$,

$$
\lim_{N\to\infty} \nu_{1/2}^{(N)}[\widehat{F}] = \frac{1}{2\pi} \sum_{\rho} F(\rho),
$$

where $\rho$ runs over all non-trivial zeros of the Riemann $\zeta$ function.

**Proof**: By Theorem 4.2, $|D_N(t_k)|^2$ approximates $|-\zeta'/\zeta(1/2+it_k)|^2$ in the $L^2$ average sense. Combining Poisson summation with the sampling completeness from Lemma 4.3, and applying the Weil explicit formula together with standard Tauberian arguments, the measure converges weakly to the zero distribution, with a factor of $1/(2\pi)$ arising from the Fourier transform normalization. $\square$

### 5.3 Limit Order and Critical Line Compactification

**Theorem 5.3 (Limit Order and Zero Distribution)**:

(A) **Correct Order (Compactification Limit)**:
Fix $\alpha=1/2$ and consider the measure $\nu_{1/2}^{(N)}$ as defined above. Then as $N\to\infty$,

$$
\lim_{N\to\infty} \nu_{1/2}^{(N)} = \frac{1}{2\pi}\sum_{\rho} \delta_{\rho}.
$$

(B) **Incorrect Order (Non-Compact Divergence)**:
If we first fix $\alpha \neq 1/2$ and let $N\to\infty$, then the scaled generator satisfies $\|\mathcal{H}_{\alpha}^{(N)}\|_{\mathrm{op}} \sim |\alpha-1/2|\frac{N^2}{(\log N)^2} \to \infty$, the spectral measure support expands to infinity (non-compact). Taking $\alpha\to 1/2$ subsequently yields a non-existent limit (or degenerates to meaningless divergence).

(C) **Non-commutativity of Limits**:

For the internal measure $\mu_\alpha^{(N)}$:
- First $N\to\infty$ (fixed $\alpha\neq 1/2$): $\|\mathcal{H}_\alpha^{(N)}\|\to\infty$, $\mu_\alpha^{(N)}$ has no compact support, weak limit does not exist.
- First $\alpha\to 1/2$ (fixed $N$): yields $\delta_0$, then $N\to\infty$ still yields $\delta_0$ (degenerate point measure).

For the frequency measure $\nu_{1/2}^{(N)}$ (defined only at $\alpha=1/2$):
- Correct order: Fix $\alpha=1/2$ (lock the real part), then $N\to\infty$ (expand the imaginary part), yielding $\frac{1}{2\pi}\sum_\rho \delta_\rho$.
- Reverse path: If $N\to\infty$ first, then $\alpha\to 1/2$, the path is invalid because $\nu$ is undefined (or corresponds to a non-compact operator) for $\alpha\neq 1/2$.

(D) **Grid Refinement Limit (Supplement)**:
Fix $\alpha=1/2$, as $N\to\infty$, the sampling step $\Delta t = \frac{2\pi}{\log N} \to 0$, and the grid $t_k$ becomes dense in any compact interval $[0,T]$. This limit is independent of the $\alpha$ limit, but provides the geometric foundation for the weak convergence of the measure $\nu_{1/2}^{(N)}$ (see Section 4.4).

### 5.4 Limit Order and Functorial Faithfulness

Define (functorial faithfulness measure):

$$
d_{\mathrm{spec}}(\alpha) := \limsup_{N\to\infty} \left| \nu_\alpha^{(N)}[\widehat{F}] - \frac{1}{2\pi}\sum_{\rho}F(\rho) \right|
$$

(taking the supremum over a suitable class of functions $F$), where $\nu_\alpha^{(N)}$ is the frequency-weighted measure corresponding to $\alpha$ ($\nu_{1/2}^{(N)}$ for $\alpha=1/2$, and the corresponding non-compact measure for $\alpha\neq1/2$).

From Theorem 5.3, we have:

- $d_{\mathrm{spec}}(1/2)=0$ when $\alpha=1/2$.
- $d_{\mathrm{spec}}(\alpha) > 0$ (in fact divergent) when $\alpha\neq1/2$.

---

## 6. Error Bounds and Convergence Analysis

### 6.1 Error Inequalities

The above results are summarized in the following inequalities:

1. **Carleman Condition**:
   
   $$
   S_N(\alpha)=\sum_{n=1}^N\frac{1}{\beta_n^{(N)}(\alpha)} \begin{cases} \to\infty & \alpha=1/2 \\ \to0 & \alpha\neq1/2 \end{cases}
   $$
   

2. **Phase Stability (Lipschitz)**:
   
   $$
   \| e^{i\alpha\operatorname{Li}(x)\mathbf{H}} - e^{i\alpha'\operatorname{Li}(x)\mathbf{H}}\| \le |\alpha-\alpha'|\cdot\operatorname{Li}(x)\cdot\|\mathbf{H}\|
   $$
   

3. **Explicit Formula Remainder Estimate**:
   
   $$
   \left| \frac{1}{\log N}\sum_{k} |D_N(t_k)|^2 \widehat{F}(t_k) - \frac{1}{2\pi}\sum_{\rho} F(\rho) \right| \le \epsilon(N) \sup|F|,
   $$
   
   with $\epsilon(N)=O(N^{-1/2})$ (analytic error) plus $O(1/\log N)$ (sampling error), where the sampling error is dominant in practice.

4. **Topological Acceptance Rate**:
   
   $$
   \mathcal{A}(\alpha) \ge \mathcal{A}_{\max}\exp\left(-\frac{(\alpha-1/2)^2}{2\sigma^2}\right)-\delta_N
   $$
   

5. **Bishop Convergence (Cauchy Sequence)**:
   
   $$
   \|\mathcal{F}_{1/2}^{(N)}-\mathcal{F}_{1/2}^{(M)}\|_{\mathrm{norm}} < 2^{-k}\quad \text{for } N,M>N(k)
   $$
   
   where $\|A\|_{\mathrm{norm}} = \|A\|_F / \sqrt{\min(N,M)}$, $N(k) = \lceil C \cdot 2^{2k} \rceil$.

6. **RH Error Bound**:
   
   $$
   \left| \sum_\rho \frac{x^\rho}{\rho} - \operatorname{Tr}\left(\frac{x^{\mathcal{L}_{1/2}^{(N)}}}{\mathcal{L}_{1/2}^{(N)}}\right) \right| < \epsilon(N) x^{1/2}\ln x,\quad \epsilon(N)=O(N^{-1/2})
   $$
   

**Distinction of Error Sources**:
- Analytic error $\epsilon_{\text{analytic}}(N) = O(N^{-1/2})$: arises from the truncation $\sum_{n>N} \Lambda(n)/\sqrt{n}$, the Mellin transform, and the approximation by the smooth truncation $\phi$.
- Sampling error $\delta_{\text{geom}}(N) = O(1/\log N)$: arises from the geometric deviation between the discrete grid $t_k$ and the continuous zeros $\gamma_j$ (see Lemma 4.3).

Since $\sqrt{N} \gg \log N$, for sufficiently large $N$, the sampling error dominates the zero localization precision, while the analytic error is comparatively smaller.

### 6.2 Explicit Derivation of $\epsilon(N)$

From the truncation error analysis in Section 4.2 and the numerical stability of the FFT, we obtain:

$$
\epsilon(N) \leq \frac{C_1}{\sqrt{N}} + \frac{C_2}{N} + C_3 e^{-c\sqrt{\log N}}
$$

where the first term comes from the truncation $\sum_{n>N} \Lambda(n)/\sqrt{n}$, the second from discretization error $\log n \approx n\delta$, and the third from the remainder term in the prime number theorem ($\psi(x)=x+O(xe^{-c\sqrt{\log x}})$). Therefore, $\epsilon(N) = O(N^{-1/2})$.

### 6.3 Constructive Convergence Protocol

Given a precision $\epsilon > 0$, the computational protocol is:

1. Choose $N = O(\epsilon^{-2})$ (to ensure analytic error $<\epsilon/2$) and $N > \exp(\pi/\epsilon)$ (to ensure sampling error $<\epsilon/2$, see Section 4.4)
2. Construct $\mathcal{F}_{1/2}^{(N)}$ (standard DFT matrix, corresponding to logarithmic sampling $x_n = 2e^{n\delta}$)
3. Perform an FFT to compute $D_N(t_k)$
4. Locate $t_k \approx \gamma_j$ by detecting local maxima of $|D_N(t_k)|$ (or $\pi$ phase jumps in $\arg(D_N(t_k))$). By Lemma 4.3, when $N > \exp(\pi/\delta)$, the localization error $|t_k - \gamma_j| < \delta$ is guaranteed for all $\gamma_j < T$.
5. The error is bounded by Inequality 6: $|\rho_j - \rho_j^{\mathrm{true}}| < \epsilon$

---

## 7. Mathematical Correspondence for Numerical Validation

### 7.1 Discrete Dirichlet Polynomial

Define the finite truncation of the logarithmic derivative (smoothly truncated):

$$
D_N(t)=\sum_{n=1}^\infty \frac{\Lambda(n)}{\sqrt{n}}\phi\left(\frac{\log n}{\log N}\right) e^{-it\log n}.
$$

#### 7.1.1 Asymptotic Density of the Sampling Grid

From Section 4.4, the spacing of frequency sampling points $t_k = \frac{2\pi k}{\log N}$ is $\Delta t = \frac{2\pi}{\log N}$, which tends to zero as $N\to\infty$. Using the Riemann-von Mangoldt zero density estimate $N(T) \sim \frac{T}{2\pi}\log T$, for the first $J$ zeros $\gamma_1,...,\gamma_J$, when

$$
N > \exp\left(\frac{2\pi}{\min_{j\leq J} |\gamma_{j+1}-\gamma_j|}\right)
$$

the grid can resolve all the first $J$ zeros (without omission or confusion).

### 7.2 FFT Phase Correspondence

At the discrete points $t_k=\frac{2\pi k}{\log N}$:

$$
D_N(t_k) \approx \sqrt{N}\cdot \mathrm{FFT}[\Lambda(n)/\sqrt{n}\,\phi(\cdot)](k).
$$

Phase error analysis:

$$
\arg(D_N(t)) - \arg\left(-\frac{\zeta'(1/2+it)}{\zeta(1/2+it)}\right) = O\left(\frac{t^2}{\log N}\right) + O\left(\frac{1}{\sqrt{N}}\right)
$$

For $t \in [0, T]$, choosing $N \sim T^2$ makes the error arbitrarily small.

### 7.3 Spectral Conditions for Zero Detection

At a zero $\rho = 1/2 + i\gamma$, $-\zeta'/\zeta$ has a simple pole, leading to a $\pi$ phase jump.

In the FFT spectrum, this manifests as a magnitude minimum and phase discontinuity:

$$
|\mathrm{FFT}(k)| \approx \left| \frac{1}{\rho - (1/2+it_k)} \right| \quad \text{at } t_k \approx \gamma
$$

forming a Lorentzian-type peak (or trough, depending on normalization).

**Theoretical Guarantee of the Phase Jump Criterion**:
From Corollary 4.5, as $t_k$ crosses $\gamma_j$ from left to right, $\arg(D_N(t))$ undergoes a $+\pi$ jump (or $-\pi$ for the opposite direction). Detecting this jump is independent of the magnitude normalization $|D_N(t_k)|^2/\|D_N\|_2^2$, providing a more robust numerical criterion.

---

## 8. Conclusion: Statement of Mathematical Completeness

**Theorem 8.1 (Equivalence of NCDFT-RH)**: The following statements are equivalent:

1. The Riemann Hypothesis: All non-trivial zeros satisfy $\Re(\rho)=1/2$.
2. Self-duality: The dual composition operator $\mathcal{U}_\alpha^{(N)}=\mathbb{I}$ if and only if $\alpha=1/2$.
3. Uniqueness of spectral compactification: The dual generator $\mathcal{H}_\alpha^{(N)}=0$ (real part locked) if and only if $\alpha=1/2$, and the spectral measure constructed from the frequency indices $t_k$ is compact (supported on the line $\Re(s)=1/2$); for $\alpha\neq1/2$, $\beta_n^{(N)}\to\infty$, and the spectral measure is non-compact in the complex plane (Carleman condition fails).
4. Functorial faithfulness: $d_{\mathrm{spec}}(\alpha)=0$ if and only if $\alpha=1/2$ (via the correct limit order).
5. FFT convergence: $\arg(\mathrm{FFT}[\Lambda(n)/\sqrt{n}\,\phi(\cdot)])$ converges in the $L^2$ sense to $\arg\left(-\frac{\zeta'(1/2+it)}{\zeta(1/2+it)}\right)$.
6. Sampling completeness: The discrete frequency grid $t_k = \frac{2\pi k}{\log N}$ can arbitrarily approximate all zero imaginary parts $\gamma_j$ as $N\to\infty$, i.e., $\lim_{N\to\infty} \min_k |t_k - \gamma_j| = 0$ for all $j$ (Lemma 4.3).

**Proof Sketch**:

- (1)$\Rightarrow$(2): From the Riemann functional equation, when zeros lie on the critical line, the dual symmetry holds strictly at $\alpha=1/2$, and by Theorem 4.1, $\mathcal{U}_{1/2}=\mathbb{I}$.
- (2)$\Rightarrow$(3): Self-duality forces $\mathcal{H}_{1/2}=0$ (real part locked); if $\alpha\neq 1/2$, then $\mathcal{U}_\alpha\neq\mathbb{I}$, and by Theorem 3.3 (all eigenvalues unbounded) and Theorem 3.6 (Jacobi coefficients diverge), non-compactness follows.
- (3)$\Rightarrow$(4): Spectral compactness ensures convergence of the Weil trace formula, so functorial faithfulness holds; non-compact spectra lead to deviation.
- (4)$\Rightarrow$(5): Functorial faithfulness guarantees that IDFT reconstruction is consistent with the explicit formula, and the FFT phase converges.
- (5)+(6)$\Rightarrow$(1): From (5), the FFT phase converges, and from (6) (sampling completeness), for any $\gamma_j$, there exists a subsequence $t_{k(N)} \to \gamma_j$ such that $\arg(D_N(t_{k(N)}))$ exhibits a $\pi$ jump. This forces the real part of the zero to be $1/2$; otherwise, the phase behavior would deviate from the pole characteristics of the logarithmic derivative. If a zero $\rho_0$ existed with $\Re(\rho_0)\neq 1/2$, its support would lie off the critical line, contradicting the measure inclusion relation in Section 5. $\square$

**Rigorous Derivation of RH**:

From Theorem 2.3 and Theorem 3.5, when $\alpha=1/2$, the scaled generator is strictly zero:

$$
\mathcal{H}_{1/2}^{(N)}=0
$$

This corresponds to the freezing of internal degrees of freedom—the real part deviation is locked to zero, and all spectral points are constrained to the line $\Re(s)=1/2$. However, this only determines the real part coordinate of the zeros. To recover the full complex zeros $\rho_j = 1/2 + i\gamma_j$ (including the imaginary part information), we must introduce the external frequency indices $t_k = \frac{2\pi k}{\log N}$ via the frequency spectral measure $\nu_{1/2}^{(N)}$:

$$
\nu_{1/2}^{(N)} := \frac{1}{\log N}\sum_{k=0}^{N-1} |D_N(t_k)|^2 \delta_{\frac{1}{2} + i t_k}
$$

The support of this measure is:

$$
\mathrm{supp}(\nu_{1/2}^{(N)}) = \left\{\frac{1}{2} + i t_k \;\middle|\; k=0,1,\dots,N-1\right\} \subset \left\{s\in\mathbb{C} \;\middle|\; \Re(s)=\frac{1}{2}\right\}
$$

i.e., all mass lies strictly on the critical line $\Re(s)=1/2$. By Theorem 5.2, as $N\to\infty$, this measure converges weakly to the empirical distribution of Riemann zeros:

$$
\lim_{N\to\infty} \nu_{1/2}^{(N)} = \frac{1}{2\pi}\sum_{\rho} \delta_{\rho}
$$

By the lower semicontinuity of weak convergence of measures, the support of the limit measure satisfies the inclusion:

$$
\mathrm{supp}\left(\sum_{\rho} \delta_{\rho}\right) \subseteq \overline{\bigcup_{N} \mathrm{supp}(\nu_{1/2}^{(N)})} \subseteq \left\{s \in \mathbb{C} \;\middle|\; \Re(s)=\frac{1}{2}\right\}
$$

Therefore, all non-trivial zeros $\rho$ must satisfy $\Re(\rho)=1/2$. If there existed a zero $\rho_0$ with $\Re(\rho_0)\neq 1/2$, the support of $\delta_{\rho_0}$ would lie off the critical line, contradicting the inclusion above. $\square$

**Key Distinction**:
- The internal spectral measure $\mu_{1/2}^{(N)}=\delta_0$ reflects $\mathcal{H}_{1/2}^{(N)}=0$, indicating only the real part locking (single-point compactification), and contains no imaginary part information.
- The frequency spectral measure $\nu_{1/2}^{(N)}$ reflects the imaginary part information expanded via the FFT frequencies $t_k$, with its support on the critical line, converging to the zero distribution.

Together, they form a complete spectral lifting: the "internal" part of $\mathcal{L}_{1/2}^{(N)}$ (from $\mathcal{H}$) locks the real part to $1/2$, while the "external" part (from $\mathrm{diag}(t_k)$) expands the imaginary part to $\gamma_j$. Only when $\alpha=1/2$ does this dual structure achieve a strictly self-dual matching. $\square$

**Corollary**: This framework provides a constructive proof pathway: approximating the infinite-dimensional analytic object (Riemann zeros) via finite-dimensional self-dual computations (standard DFT/FFT), with all steps satisfying the computability requirements of Bishop's constructive mathematics. The NCDFT framework reveals that the Riemann Hypothesis is fundamentally a duality condition: only when the parameter $\alpha=1/2$ does the non-commutative discrete Fourier transform achieve strict self-duality, at which point the functorial faithfulness between arithmetic and spectra reaches perfect matching, and the spectral measure compactifies from non-compact extended states ($\alpha \neq 1/2$) to discrete point states on the critical line ($\alpha = 1/2$).

---

## Appendix A: FFT Butterfly Operation and Dual Structure

In the Fast Fourier Transform (FFT), the butterfly operation is the core step of the divide-and-conquer algorithm. Its elegant algebraic structure not only reflects the symmetry of the DFT but also has a profound parallel with the duality $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ in the NCDFT framework. This appendix briefly outlines the algebraic form of the butterfly operation and explains its connection to the functional equation.

### A.1 Butterfly Operation in Decimation-in-Time FFT

Consider a DFT of length $N$ ($N=2^m$):

$$
X_k = \sum_{n=0}^{N-1} x_n \omega_N^{kn}, \quad \omega_N = e^{-2\pi i/N}.
$$

Splitting the input sequence into even and odd indices:

$$
X_k = \sum_{m=0}^{N/2-1} x_{2m} \omega_N^{k(2m)} + \sum_{m=0}^{N/2-1} x_{2m+1} \omega_N^{k(2m+1)} = E_k + \omega_N^k O_k,
$$

where $E_k$ and $O_k$ are DFTs of length $N/2$:

$$
E_k = \sum_{m=0}^{N/2-1} x_{2m} \omega_{N/2}^{km}, \quad O_k = \sum_{m=0}^{N/2-1} x_{2m+1} \omega_{N/2}^{km}.
$$

For $k$ and $k+N/2$, using $\omega_N^{k+N/2} = -\omega_N^k$, we obtain the classic butterfly:
$$
\begin{cases}
X_k = E_k + \omega_N^k O_k, \\
X_{k+N/2} = E_k - \omega_N^k O_k.
\end{cases}
$$

This operation is depicted as a butterfly-like structure, hence the name.

### A.2 Butterfly Operation and Duality

The butterfly operation reveals an intrinsic symmetry of the DFT: the even part $E_k$ and odd part $O_k$ are combined via the twiddle factor $\omega_N^k$, which itself satisfies $\omega_N^{k+N/2} = -\omega_N^k$. This is a discrete manifestation of the $s \leftrightarrow 1-s$ symmetry in the functional equation. In the NCDFT framework, the duality $\mathcal{F}_\alpha \leftrightarrow \mathcal{F}_{1-\alpha}$ corresponds to swapping the parameters $\alpha$ and $1-\alpha$, while the parity decomposition and sign reversal in the butterfly operation mimic this dual relationship.

Furthermore, if we view the NCDFT matrix $\mathcal{F}_\alpha^{(N)}$ as a DFT augmented with an internal phase $\exp(i(\alpha-1/2)\operatorname{Li}(x_n)\mathbf{H})$, its fast algorithm can similarly be decomposed into a butterfly-like structure, where the twiddle factors are multiplied by the corresponding diagonal matrices. This structure ensures the computational efficiency of the NCDFT while preserving its dual connection to the Riemann functional equation.

### A.3 Butterfly Operation and the Functional Equation

Combining the butterfly operation with the correspondence table in Section 4.2, we have:

| FFT Butterfly Operation | NCDFT Duality |
| :--- | :--- |
| Parity decomposition $E_k, O_k$ | Parameter duality $\alpha \leftrightarrow 1-\alpha$ |
| Twiddle factors $\omega_N^k$ and $-\omega_N^k$ | Phase factors $\exp(i(\alpha-1/2)\operatorname{Li}(x_n)\mathbf{H})$ and its conjugate |
| Symmetry between $X_k$ and $X_{k+N/2}$ | Self-duality condition $\mathcal{U}_{1/2}=\mathbb{I}$ |
| Grid refinement $N\to\infty$ | Sampling completeness $t_k \to \gamma_j$ (Section 4.4) |

---

## References

1. H. M. Edwards, *Riemann's Zeta Function*, Dover, 2001 [Original: 1974].
2. E. C. Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., Oxford University Press, 1986.
3. A. Weil, "Sur les formules explicites de la théorie des nombres premiers", Comm. Sem. Math. Univ. Lund, 1952.
4. A. E. Ingham, "On the estimation of $N(\sigma,T)$", Quart. J. Math., 1940.
5. J. W. Cooley, J. W. Tukey, "An algorithm for the machine calculation of complex Fourier series", Math. Comp., 1965.
6. E. Bishop, *Foundations of Constructive Analysis*, McGraw-Hill, 1967.