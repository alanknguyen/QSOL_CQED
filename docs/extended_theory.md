# Extended theory notes

Supplementary derivations that go beyond the manuscript ([../paper/qsol_cqed.pdf](../paper/qsol_cqed.pdf)): the Husimi $Q$-function, the cat-state Wigner function and fidelity, the enhanced decoherence rate, photon blockade, Bloch-vector geometry, the photon-number parity comb, and the interpretation of the cat-survival phase diagram. Numbers quoted here are those of the corrected simulations.


> **Note:** The full theoretical framework — JC Hamiltonian derivation, rotating-wave approximation, dressed states and eigenvalues, collapse and revival timescales, Wigner function formalism, Lindblad master equation, and von Neumann entropy — is presented in the [companion paper](../paper/qsol_cqed.pdf) (Sections II–III). The sections below cover **new theoretical material** developed for the extended computational results in this repository.

### 1. Husimi Q-Function and the Phase-Space Smoothing Theorem

The Husimi Q-function is defined as the diagonal matrix element of the density operator in the coherent-state basis:

$$Q(\alpha) = \frac{1}{\pi}\langle\alpha|\rho|\alpha\rangle$$

Since $\rho$ is a positive operator and $|\alpha\rangle$ is a normalized state, it follows immediately that $Q(\alpha) \geq 0$ for all $\alpha$. This is in sharp contrast to the Wigner function, which can take negative values.

The precise relationship between $Q$ and $W$ is given by the **Gaussian convolution theorem**. In the phase-space convention used throughout this repository, $a = (x + ip)/\sqrt{2}$, a coherent state $|\alpha\rangle$ is centred at $(x, p) = \sqrt{2}(\mathrm{Re}\,\alpha, \mathrm{Im}\,\alpha)$ and the vacuum Wigner function is $W_{|0\rangle}(x,p) = \frac{1}{\pi}e^{-(x^2 + p^2)}$ (variance $1/2$ per quadrature). Then

$$Q(x, p) = \frac{1}{\pi}\int\!\!\int W(x', p')\,\exp\!\left[-\left((x' - x)^2 + (p' - p)^2\right)\right]dx'\,dp' = (W * W_{|0\rangle})(x,p),$$

normalized so that $\int Q\,dx\,dp = 1$. In other words, $Q$ is obtained by **smoothing $W$ with the vacuum Wigner function**, a minimum-uncertainty Gaussian of width equal to the vacuum fluctuation.

This smoothing has an irreversible information-theoretic consequence: a feature of $W$ with wavenumber $k$ is attenuated by the Fourier transform of the kernel, $e^{-k^2/4}$. For a cat state with components at $\pm\beta$ (separation $2\sqrt{2}|\beta|$ in these units) the interference term of $W$ oscillates with wavenumber $k = 2\sqrt{2}|\beta|$ (fringe spacing $\pi/(\sqrt{2}|\beta|)$; for $\bar{n} = |\beta|^2 = 10$ this is $0.70$, as observed). Convolving $e^{-p^2}\cos(kp)$ with the kernel gives an oscillation of *half* the wavenumber with amplitude reduced by $e^{-k^2/8} = e^{-|\beta|^2} = e^{-\bar{n}}$. For $\bar{n} = 10$ the residual fringes in $Q$ are suppressed by $e^{-10} \approx 5\times10^{-5}$, which is why the Q-function shows two smooth blobs while the Wigner function resolves the oscillatory structure (as demonstrated in our Wigner vs Q comparison figures). The same factor follows directly from $Q(\beta') \propto |\langle\beta'|\text{cat}\rangle|^2$: the cross term at the midpoint $\beta' = 0$ is $e^{-|\beta|^2}$ relative to the lobe maxima.

This can also be understood information-theoretically through the **Wehrl entropy** $h_W = -\int Q(\alpha)\ln Q(\alpha)\,d^2\alpha$, which satisfies the Lieb bound $h_W \geq 1$ (with equality only for coherent states) and is never smaller than the von Neumann entropy, $h_W \geq S(\rho)$; the gap quantifies the information lost by the Gaussian smoothing.

The practical implication is that **Wigner negativity, not Q-function structure, is the correct witness of quantum coherence in phase space.** The Q-function is useful for visualization but fundamentally cannot distinguish a quantum cat state $|\alpha\rangle + |-\alpha\rangle$ from a classical mixture $|\alpha\rangle\langle\alpha| + |-\alpha\rangle\langle-\alpha|$.

### 2. Schrödinger Cat-State Wigner Function and Fidelity

At time $t = t_r/2$, the JC interaction produces a cavity field state that closely approximates a two-component cat state whose branches sit at $\beta = \pm i\sqrt{\bar{n}}$ (on the $p$ axis, rotated by $\pm\pi/2$ from the initial amplitude $\alpha = \sqrt{\bar{n}}$) with a relative phase $\theta \approx 1.2\pi$ for $\bar{n} = 10$ — closer to an *odd* than to an even cat (photon-number parity $\langle(-1)^{\hat n}\rangle = -0.64$). For the algebra below we take the even cat along the real axis as the reference example; the JC cat is obtained by rotating $x \leftrightarrow p$ and changing the relative phase. The ideal even cat state and its normalization are

$$|\text{cat}^+\rangle = \mathcal{N}_+\bigl(|\alpha\rangle + |-\alpha\rangle\bigr), \qquad \mathcal{N}_+ = \frac{1}{\sqrt{2(1 + e^{-2|\alpha|^2})}}$$

The Wigner function of this state decomposes exactly into three terms (convention $a = (x+ip)/\sqrt{2}$):

$$W_\text{cat}(x,p) = \mathcal{N}_+^2\Bigl[W_\alpha(x,p) + W_{-\alpha}(x,p)\Bigr] + 2\mathcal{N}_+^2\,W_\text{int}(x,p)$$

where $W_{\pm\alpha}$ are the Gaussian Wigner functions of the two coherent components, centred at $x = \pm\sqrt{2}\alpha$, and the **interference term** is

$$W_\text{int}(x,p) = \frac{1}{\pi}\cos\!\bigl(2\sqrt{2}\,\alpha\,p\bigr)\,e^{-(x^2 + p^2)}$$

for $\alpha$ real ($\int W_\text{int}\,dx\,dp = e^{-2\alpha^2}$, which fixes the normalization). This interference term is responsible for the oscillatory fringes between the two lobes and is the sole source of Wigner negativity. The fringe spacing along the axis perpendicular to the branch separation is

$$\Delta_\text{fringe} = \frac{\pi}{\sqrt{2}\,\alpha} = \frac{\pi}{\sqrt{2\bar{n}}}$$

which for $\bar{n} = 10$ gives $0.70$, consistent with the simulated cross-section $W(x, p{=}0)$ (there the branches lie on the $p$ axis, so the fringes run along $x$).

The most negative value of the ideal cat is reached at the first dark fringe next to the origin, $2\sqrt{2}\alpha p = \pi$, where $W_\text{cat} \approx -\tfrac{2}{\pi}e^{-\pi^2/(8\alpha^2)} \approx -0.28$ for $\bar{n} = 10$. The simulated JC cat, which has fidelity 0.78 with the best two-component cat, reaches $W \approx -0.25$.

A key result is that **larger cat states have finer fringes.** Since $\Delta p \propto 1/\sqrt{\bar{n}}$, the fringe frequency grows with $\bar{n}$, making the interference pattern increasingly vulnerable to any smoothing process — whether instrumental (finite detector resolution) or physical (decoherence). This is the phase-space manifestation of the fragility of macroscopic superpositions.

We quantify the quality of the dynamically generated cat state using the **fidelity** with the best-fitting two-component cat,

$$F_\text{cat}(\bar{n}, \kappa/g) = \max_{\beta,\,\theta}\;\langle\text{cat}(\beta,\theta)|\,\rho_\text{field}(t_r/2)\,|\text{cat}(\beta,\theta)\rangle, \qquad |\text{cat}(\beta,\theta)\rangle \propto |\beta\rangle + e^{i\theta}|-\beta\rangle,$$

with $\beta$ scanned in a neighbourhood of $i\sqrt{\bar{n}}$ and $\theta \in [0, 2\pi)$. Comparing instead with the even cat $|\alpha\rangle + |-\alpha\rangle$ on the real axis gives $F \approx 0$ for every $\bar{n} \gtrsim 3$, simply because the JC branches sit on the $p$ axis with a non-even relative phase — a comparison with the wrong reference state, not a property of the dynamics. With the correct reference, our phase diagram shows $F_\text{cat} \approx 0.78$ at $\kappa = 0$ for every $\bar{n}$ from 2 to 25: the anharmonic $\sqrt{n+1}$ spectrum limits the cat quality to a fixed value rather than degrading it with $\bar{n}$. Dissipation lowers $F_\text{cat}$ (0.54 and 0.46 at $\kappa/g = 0.01$ and 0.02 for $\bar{n} = 10$), and the $F_\text{cat} = 0.5$ boundary scales as $\bar{n}^{-1.25}$, tracking the $\delta = 0.05$ boundary.

### 3. Microscopic Derivation of Enhanced Decoherence

The paper quotes the fringe decoherence rate $\Gamma_\text{dec} = 2\kappa\bar{n}$ following Zurek [22]. Here we derive it explicitly from the Lindblad dynamics.

Consider the master equation $\dot{\rho} = \kappa\mathcal{D}[a]\rho$ acting on the cavity field alone (ignoring the atom for this argument). We wish to compute the decay rate of the off-diagonal coherence $\langle\alpha|\rho|-\alpha\rangle$ between two coherent states separated by $2\alpha$ in phase space.

Using the coherent-state matrix elements of the Lindblad superoperator, and the eigenvalue relation $a|\alpha\rangle = \alpha|\alpha\rangle$:

$$\frac{d}{dt}\langle\alpha|\rho|-\alpha\rangle = \kappa\left[\alpha(-\alpha)^* - \frac{|\alpha|^2 + |-\alpha|^2}{2}\right]\langle\alpha|\rho|-\alpha\rangle$$

$$= \kappa\left[-|\alpha|^2 - |\alpha|^2\right]\langle\alpha|\rho|-\alpha\rangle = -2\kappa|\alpha|^2\,\langle\alpha|\rho|-\alpha\rangle$$

This gives the exact result:

$$\langle\alpha|\rho(t)|-\alpha\rangle = \langle\alpha|\rho(0)|-\alpha\rangle\,\exp\!\left(-2\kappa|\alpha|^2\,t\right)$$

The decoherence rate is therefore

$$\boxed{\Gamma_\text{dec} = 2\kappa|\alpha|^2 = 2\kappa\bar{n}}$$

This is a factor of $2\bar{n}$ faster than the energy decay rate $\kappa$. The physical interpretation is that each photon lost from the cavity carries "which-component" information about the cat state: a measurement of the photon's phase would distinguish $|\alpha\rangle$ from $|-\alpha\rangle$. The rate of such information leakage is proportional to the photon number.

The **decoherence time** is

$$t_\text{dec} = \frac{1}{2\kappa\bar{n}}$$

For the cat state to be observable, we need $t_\text{dec} \gtrsim t_r/2$, i.e., the cat must survive long enough to form. This requires

$$\frac{1}{2\kappa\bar{n}} \gtrsim \frac{\pi\sqrt{\bar{n}}}{g} \qquad \Longrightarrow \qquad \frac{\kappa}{g} \lesssim \frac{1}{2\pi\bar{n}^{3/2}}$$

For $\bar{n} = 10$ this crude criterion gives $\kappa/g \lesssim 0.005$, whereas the simulated $\delta = 0.05$ contour lies at $\kappa/g \approx 0.034$: the estimate is conservative because the branch separation, and with it the decoherence rate, builds up gradually during the collapse rather than being $2\sqrt{\bar{n}}$ from $t = 0$. The measured boundary scales as $\kappa_c/g \approx 0.46\,\bar{n}^{-1.14}$ (fit over $\bar{n} = 4$–$25$; $\bar{n}^{-1.26}$ for $\bar{n} \geq 8$, and $\bar{n}^{-1.25}$ for the $F_\text{cat} = 0.5$ boundary), i.e. between the naive $1/\bar{n}$ and the $\bar{n}^{-3/2}$ of this argument, which combines the $\bar{n}$-dependent decoherence rate with the $\bar{n}^{1/2}$-dependent formation time.

### 4. Second-Order Coherence and the Photon Blockade Effect

The equal-time second-order correlation function is defined as

$$g^{(2)}(0) = \frac{\langle a^\dagger a^\dagger a\,a\rangle}{\langle a^\dagger a\rangle^2} = \frac{\langle \hat{n}(\hat{n}-1)\rangle}{\langle \hat{n}\rangle^2} = 1 + \frac{\mathrm{Var}(\hat{n}) - \langle\hat{n}\rangle}{\langle\hat{n}\rangle^2}$$

This classifies photon statistics:

| $g^{(2)}(0)$ | Statistics | Physical meaning |
|:---:|:---:|:---|
| $0$ | Perfect antibunching | At most one photon present (single-photon source) |
| $< 1$ | Sub-Poissonian | Photon-number variance below shot noise |
| $= 1$ | Poissonian | Coherent-state (laser) statistics |
| $= 2$ | Thermal (chaotic) | Bose-Einstein bunching |
| $> 2$ | Super-thermal | Exotic multi-photon correlations |

The connection to the paper's Sec. V is that sub-Poissonian statistics ($g^{(2)}(0) < 1$) and photon antibunching are among the definitive criteria for non-classical light. Here we develop the full theory of how the JC nonlinearity produces antibunching via the photon blockade.

#### The Driven Dissipative JC Model

Adding a coherent drive of amplitude $\varepsilon$ at frequency $\omega_L$ and moving to the rotating frame:

$$H_\text{driven} = \Delta_c\,a^\dagger a + \frac{\Delta_a}{2}\sigma_z + g\bigl(a^\dagger\sigma^- + a\,\sigma^+\bigr) + \varepsilon\bigl(a^\dagger + a\bigr)$$

where $\Delta_c = \omega_c - \omega_L$ and $\Delta_a = \omega_a - \omega_L$. The steady state $\rho_\text{ss}$ is found from the Lindblad equation $0 = -i[H, \rho_\text{ss}] + \kappa\mathcal{D}[a]\rho_\text{ss} + \gamma\mathcal{D}[\sigma^-]\rho_\text{ss}$, and $g^{(2)}(0)$ is computed directly from $\rho_\text{ss}$.

#### The Blockade Mechanism

The photon blockade arises from the **anharmonicity of the JC energy ladder**. On resonance ($\Delta_c = \Delta_a = 0$), the dressed-state energies (from the paper's Eq. 14) give transition frequencies:

$$\omega_{0 \to 1,\pm} = \omega_c \pm g \qquad \text{(first photon)}$$

$$\omega_{1,\pm \to 2,\pm} = \omega_c \pm g(\sqrt{2} - 1) \quad \text{(same branch)}, \qquad \omega_{1,\mp \to 2,\pm} = \omega_c \pm g(\sqrt{2} + 1) \quad \text{(cross branch)} \qquad \text{(second photon)}$$

The energy cost of the first photon differs from the second by

$$\delta E_\text{block} = \hbar g\bigl(2 - \sqrt{2}\bigr) \approx 0.59\,\hbar g$$

A drive laser tuned to the first transition (say $\omega_L = \omega_c - g$) is therefore **detuned from the second transition** by $\delta E_\text{block}/\hbar$. If this detuning exceeds the dressed-state linewidth ($\sim \kappa$), the second photon cannot be absorbed while the first is still in the cavity. This is the photon blockade.

The blockade condition is:

$$\frac{g(2 - \sqrt{2})}{\kappa} \gg 1 \qquad \Longleftrightarrow \qquad \frac{g}{\kappa} \gg \frac{1}{2 - \sqrt{2}} \approx 1.71$$

In practice, our simulations (drive on the polariton, $\varepsilon/\kappa = 0.05$) give $g^{(2)}(0) = 0.11$ at $g/\kappa = 5$ and $0.03$ at $g/\kappa = 10$, the finite linewidth and multi-level structure softening the ideal step. The laser must sit on a polariton for this to work: on bare-cavity resonance ($\omega_L = \omega_c$) it is detuned by $g$ from both single-photon transitions but only by $g/\sqrt{2}$ from the two-photon transition $|0\rangle \to |2,\pm\rangle$, so the cavity is nearly empty ($\langle n\rangle \sim 10^{-8}$ at $g/\kappa = 5$) and $g^{(2)}(0) \sim 10^7$ — bunching, not blockade.

#### Drive Strength Crossover

At fixed $g/\kappa$, increasing the drive $\varepsilon$ eventually overwhelms the blockade. In the limit $\varepsilon \gg g$, the cavity is driven classically and $g^{(2)}(0) \to 1$. In our sweeps ($\Delta = g$, $\varepsilon/\kappa \leq 3$) the rise is monotonic: $g^{(2)}(0)$ approaches 1 from below, more slowly for larger $g/\kappa$.

#### Photon Blockade Spectrum

Sweeping the laser detuning $\Delta$ at fixed $g, \kappa, \varepsilon$ maps out the spectral structure of the blockade. The transmission spectrum $\langle n\rangle(\Delta)$ shows the vacuum Rabi doublet (peaks at $\Delta = \pm g$), while $g^{(2)}(0)(\Delta)$ reveals:

- **Antibunching dips** *at* the polariton resonances $\Delta = \pm g$ ($g^{(2)}(0) \approx 0.1$ for $g/\kappa = 5$), where the first photon is absorbed resonantly and the second is blocked
- **Bunching spikes** at the two-photon resonances $\Delta = \pm g/\sqrt{2}$, where $2\omega_L$ matches $E_{2,\pm} = 2\omega_c \pm \sqrt{2}g$ ($g^{(2)}(0) \sim 10^2$), and a very large bunching peak on bare resonance $\Delta = 0$, where the single-photon transitions are far off resonance and the cavity is almost empty

This spectral structure is the fingerprint of the quantized JC ladder and has been directly observed in circuit QED experiments.

### 5. Bloch Vector Dynamics and the Geometry of Entanglement

The reduced atomic state is a $2 \times 2$ density matrix, completely characterized by the Bloch vector:

$$\mathbf{r} = \bigl(\langle\sigma_x\rangle,\, \langle\sigma_y\rangle,\, \langle\sigma_z\rangle\bigr) = \bigl(\mathrm{Tr}[\rho_\text{atom}\sigma_x],\, \mathrm{Tr}[\rho_\text{atom}\sigma_y],\, \mathrm{Tr}[\rho_\text{atom}\sigma_z]\bigr)$$

A pure atomic state lies on the surface of the Bloch sphere ($|\mathbf{r}| = 1$); a maximally mixed state sits at the center ($|\mathbf{r}| = 0$). The Bloch vector length is related to the atomic purity by

$$\mathrm{Tr}[\rho_\text{atom}^2] = \frac{1 + |\mathbf{r}|^2}{2}$$

and to the von Neumann entropy by the binary entropy function:

$$S(\rho_\text{atom}) = h\!\left(\frac{1 + |\mathbf{r}|}{2}\right), \qquad h(p) = -p\log_2 p - (1{-}p)\log_2(1{-}p)$$

This provides a geometric interpretation of entanglement dynamics: **the atom's entanglement with the field is encoded in how far the Bloch vector has retreated from the surface toward the center** (for a globally pure state). During the collapse $|\mathbf{r}|$ drops rapidly from 1 to $\approx 0.2$ as the atom becomes strongly entangled with the field. At the cat-state time $t = t_r/2$ the Bloch vector has returned to $|\mathbf{r}| = 0.96$, close to the surface: the atom is nearly pure again because it has **disentangled** from the field, which is left in a nearly pure cat state (Gea-Banacloche, PRL **65**, 3385 (1990)). The cat state and the entanglement minimum are two descriptions of the same event.

The trajectory itself is not a simple radial contraction. In 3D, the Bloch vector traces a **spiral** because the coherent Rabi dynamics (rotation about an axis in the $xz$-plane) compete with the dephasing caused by the spread of Fock-state Rabi frequencies. The spiral structure is visible in our Bloch sphere animations: at early times the atom precesses rapidly (Rabi oscillations), while the envelope of the spiral contracts (collapse). The spiral re-expands toward the surface at $t_r/2$, then contracts again around the revival $t_r$ ($|\mathbf{r}| \approx 0.5$, oscillating at the Rabi frequency) as the branches rephase and re-entangle.

### 6. Photon Number Distribution Dynamics and the Number-Space Cat Signature

The photon number distribution $P(n,t) = \langle n|\rho_\text{field}(t)|n\rangle$ provides complementary information to the Wigner function. For the initial state $|e\rangle \otimes |\alpha\rangle$ evolving under the resonant JC Hamiltonian, each component $|e,n\rangle$ evolves into $\cos(g\sqrt{n+1}\,t)|e,n\rangle - i\sin(g\sqrt{n+1}\,t)|g,n{+}1\rangle$, so the reduced field has the exact diagonal elements

$$P(n,t) = p(n)\cos^2\!\bigl(g\sqrt{n+1}\,t\bigr) + p(n-1)\sin^2\!\bigl(g\sqrt{n}\,t\bigr), \qquad p(n) = e^{-\bar{n}}\frac{\bar{n}^n}{n!}.$$

At $t = 0$ this is the Poisson distribution. At $t = t_r/2 = \pi\sqrt{\bar{n}}/g$ the phases $g\sqrt{n+1}\,t \approx \pi\bar{n} + \tfrac{\pi}{2}(n + 1 - \bar{n})$ alternate by $\pi/2$ from one $n$ to the next, so $\cos^2$ and $\sin^2$ switch between 0 and 1 with the *parity* of $n$: $P(n,t_r/2)$ becomes a comb in which (for $\bar{n} = 10$) odd photon numbers carry 82% of the weight and $\langle(-1)^{\hat n}\rangle = -0.64$. This is the number-space signature of the cat state: the two branches at $\pm i\sqrt{\bar{n}}$ have the *same* $|\beta|$, so their Poisson envelopes coincide and no bimodal split appears; what survives is their interference, an (approximately odd) parity comb. Such a comb provides an experimentally accessible signature of cat-state parity that does not require full Wigner tomography.

At the revival time $t_r$ the phases alternate by $\pi$ instead, the comb disappears and a single hump returns, but the anharmonicity of the JC spectrum ($\sqrt{n+1}$ rather than linear) prevents perfect recurrence. The distribution at $t = t_r$ is broader and more irregular than the initial Poissonian, with a variance that exceeds $\bar{n}$ — the field has acquired super-Poissonian statistics through its interaction with the atom.

### 7. Cat-State Survival Phase Diagram: Interpretation

The 2D parameter sweep over $(\bar{n}, \kappa/g)$ is the central new computational result. The phase diagram encodes two competing effects:

**Negativity grows with $\bar{n}$ at $\kappa = 0$.** In the ideal (lossless) case, larger $\bar{n}$ produces cat states with sharper interference fringes and larger $\delta$. The negativity saturates near $\delta \approx 0.88$ for $\bar{n} \gtrsim 14$, while the cat fidelity stays at $F_\text{cat} \approx 0.78$ for all $\bar{n} \geq 2$.

**Decoherence accelerates with $\bar{n}$ at $\kappa > 0$.** From Sec. 3, the decoherence rate $\Gamma_\text{dec} = 2\kappa\bar{n}$ grows linearly with $\bar{n}$, while the formation time $t_r/2 = \pi\sqrt{\bar{n}}/g$ grows as $\sqrt{\bar{n}}$. The product $\Gamma_\text{dec} \cdot t_r/2 = 2\pi\kappa\bar{n}^{3/2}/g$ grows as $\bar{n}^{3/2}$, meaning that larger cat states are **exponentially harder to observe** at any fixed nonzero $\kappa$.

The **critical decay rate** for cat-state survival scales as

$$\frac{\kappa_\text{crit}}{g} \sim \frac{C}{\bar{n}^{3/2}}$$

where $C$ is a threshold constant that depends on the observability criterion for $\delta$. Our phase diagram gives $\kappa_c/g \approx 0.46\,\bar{n}^{-1.14}$ over $\bar{n} = 4$–$25$ ($\bar{n}^{-1.26}$ for $\bar{n} \geq 8$): steeper than the naive $1/\bar{n}$ but shallower than this $\bar{n}^{-3/2}$ estimate, because the branch separation grows gradually during the collapse.

The parameter slices provide practical design guidance. To observe a cat state with $\delta > 0.05$ at $\bar{n} = 10$, one requires $\kappa/g \lesssim 0.034$, i.e., $g/\kappa \gtrsim 30$. This condition is achievable in microwave cavity QED ($g/\kappa \sim 10^4$) and circuit QED ($g/\kappa \sim 100\text{–}300$), but remains challenging in optical cavities where typical $g/\kappa \sim 1\text{–}10$.

---

