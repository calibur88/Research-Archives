# Diagrammatic Number Theory

## 1. Foundational Axioms

Let $\mathbb{F}$ be the base field, $S^1 = \mathbb{R}/2\pi\mathbb{Z}$.

**Axiom 1.1 (Phase Space Ontology)**  
Objects inhabit specific phase spaces, whose dimension constitutes the intrinsic measure of existence modalities.

**Axiom 1.2 (Temporal Polarization)**  
The temporal charge $\tau \in \mathbb{Z}_2 = \{+, -\}$ is absolutely conserved.

**Axiom 1.3 (Opposite-Direction Prohibition)**  
All operations between oppositely-oriented objects are forbidden, except for strict dual pairs $(a, a^\wedge)$:

$$
\forall a \in \mathcal{A}_+, \; b \in \mathcal{A}_-,\; a \neq b^\wedge \implies a \star b = \varnothing \quad (\star \in \{+,-,\times\})
$$

**Axiom 1.4 (Structural Uniqueness)**  
The phase state $\sigma \in \{\text{ch}, \text{ord}, \text{hom}\}$ is an ontological invariant; transitions occur solely through the zero-phase $\mathbf{0}$.

---

## 2. Hierarchical Construction

### Level 0: Singularity $\mathbf{0}$

- No diagram, $\dim = 0$, no symplectic structure
- Absorption: $\forall x,\; x \star \mathbf{0} = \mathbf{0}$
- Ubiquity: $\mathbf{0} \in \mathcal{A}_+ \cap \mathcal{A}_-$ (the unique element belonging simultaneously to positive and negative worlds)
- Nullification-Rebirth: $\mathcal{S}_i \to \mathbf{0} \to \mathcal{S}_j$ (no causal continuity)

### Level 1: Diagrammatic Elements $a \in \text{Sub}(\mathcal{I})$

Inhabits the diagram $\mathcal{I} = \mathbb{A}^2_{\omega}$ (2D symplectic space, $\omega = dq \wedge dp$):

- $(v, \phi, \tau) \in \mathbb{F}^\times \times S^1 \times \mathbb{Z}_2$
- Duality: $v(a^\wedge) = v(a)^{-1},\; \phi(a^\wedge) = \phi(a) + \pi/2,\; \tau(a^\wedge) = -\tau(a)$
- Symplectic pairing: $\omega(a, a^\wedge) = 1$ (automatically satisfied)

### Level 2: Spectra $\mathcal{S}_\sigma$ (Symplectic fibration $\pi: \mathcal{E} \to \mathcal{B}$)$

| Phase State | Symplectic Characteristic | Statistical Characteristic |
|-------------|---------------------------|----------------------------|
| $\mathcal{S}_{\text{ch}}$ | Non-trivial monodromy, ergodic | Exponentially random |
| $\mathcal{S}_{\text{ord}}$ | Trivial fibration, integrable | Periodically deterministic |
| $\mathcal{S}_{\text{hom}}$ | Statistical averaging | Coarse-grained equivalence |

### Level 3: Groups $\mathcal{G}_{\mathcal{S}_\sigma}$

- Spectrum binding: acts strictly on objects inhabiting $\mathcal{S}_\sigma$
- Opposite-direction prohibition: interaction forbidden if group element $g$ and object $x$ satisfy $\tau(g) \neq \tau(x)$
- No cross-spectrum property: no homomorphism exists $\mathcal{G}_{\mathcal{S}_i} \to \mathcal{G}_{\mathcal{S}_j}$ $(i \neq j)$

### Level 4: Meta-Objects $\mathcal{M}$

Inhabits $J^\infty(\Sigma_3)$, performing static classification and validity checks without executing transformations.

---

## 3. Operational Algebra

### Addition and Same-Direction Subtraction

$$
+: \mathcal{A}_\tau^\sigma \times \mathcal{A}_\tau^\sigma \to \mathcal{A}_\tau^\sigma \cup \{\mathbf{0}\}
$$

- Opposite-direction prohibition (non-dual)
- Self-subtraction nullification: $a - a = \mathbf{0}$

### Multiplication (Strict Dual Channel)

$$
\times: \mathcal{A}_+^\sigma \times \mathcal{A}_-^\sigma \to \mathbb{K} \cup \{\mathbf{0}\} \cup \{\varnothing\}
$$

- Defined only for $(a, a^\wedge)$ (other opposite-direction pairs forbidden)
- Symplectic volume: $v(a \times a^\wedge) = \omega(a, a^\wedge) = 1$
- Symplectic orthogonal annihilation: $\omega(a, b) = 0 \implies a \times b = \mathbf{0}$

### Division (Same-Direction Scaling)

$$
/: \mathcal{A}_\tau^\sigma \times (\mathcal{A}_\tau^\sigma)^\times \to \mathcal{A}_\tau^\sigma \cup \{\varnothing\}
$$

- Same-direction requirement: $\tau(b) = \tau(a)$, otherwise $\varnothing$
- Same-spectrum requirement, division by zero forbidden
- Independent definition: $v(b/a) = v(b)/v(a),\; \phi(b/a) = \phi(b) - \phi(a)$

---

## 4. Core Constraints

- **Opposite-Direction Operations:** All forbidden except strict dual multiplication.
- **Spectrum Consistency:** Operations between objects of different spectra yield $\varnothing$ (except through rebirth via $\mathbf{0}$).
- **Group Internality:** Group operations do not cross spectra; group elements and objects of opposite orientation cannot interact.
- **Zero-Phase Has No Diagram:** $\mathbf{0}$ is a pure singularity, not a symplectic manifold.
