# Figure and animation gallery

Captions for every animation and repository figure, with the physics as corrected in September 2026. Paper figures are described in the manuscript itself ([../paper/qsol_cqed.pdf](../paper/qsol_cqed.pdf)). Conventions: $\delta = \int|W|\,dx\,dp - 1$, $a = (x+ip)/\sqrt{2}$, $t_c \equiv 1/g$, $t_r = 2\pi\sqrt{\bar n}/g$.

---

## Animations

### Wigner Function Evolution

<p align="center">
  <img src="../animations/anim_wigner_evolution.gif" width="750">
</p>

A coherent state $|\alpha = \sqrt{10}\rangle$ evolves under the resonant JC Hamiltonian. During the collapse of Rabi oscillations, the intracavity field splits into a superposition of two phase-space components — a Schrödinger cat state. **Left:** Wigner function $W(x,p)$ of the reduced cavity field state $\rho_\text{field} = \mathrm{Tr}_\text{atom}[\rho]$, computed on a 200×200 phase-space grid. Interference fringes between the two coherent components produce negative regions ($W < 0$, blue), the hallmark of non-classicality. **Right:** Atomic inversion $\langle\sigma_z\rangle(t)$ with a moving time marker. The animation pauses at the cat-state time ($t = t_r/2$), where the atom and field approximately **disentangle** and the field is a nearly pure cat (purity 0.96, $\delta = 0.85$), and at the first revival ($t = t_r = 2\pi\sqrt{\bar{n}}/g$), where the branches rephase and the atom and field are strongly entangled again.

### Entanglement, Inversion, and Purity

<p align="center">
  <img src="../animations/anim_entanglement.gif" width="850">
</p>

Simultaneous evolution of three complementary observables. **Left:** Wigner function $W(x,p)$. **Center:** Atomic inversion $\langle\sigma_z\rangle(t)$. **Right:** Von Neumann entanglement entropy $S(\rho_\text{atom})$ (red) and field-state purity $\mathrm{Tr}[\rho_\text{field}^2]$ (green). Entropy rises to 1 bit during the collapse (atom and field maximally entangled) while the field purity drops to ~0.5: the reduced field is then a statistical mixture of two near-orthogonal branches. At $t = t_r/2$ the entropy falls to its minimum ($\approx 0.14$ bit) and the purity peaks at 0.96: the atom and field disentangle and the field is left in a nearly pure cat state. Toward the revival the entanglement grows again ($S \approx 0.8$ bit at $t_r$).

### Decoherence Destroys the Cat State

<p align="center">
  <img src="../animations/anim_decoherence.gif" width="750">
</p>

Cavity photon loss via the Lindblad dissipator $\kappa \mathcal{D}[a]$ erases quantum coherence on a timescale $1/(2\kappa \bar{n})$, far shorter than the bare cavity lifetime $1/\kappa$. **Left:** Wigner function at the cat-state time $t = t_r/2$ as the decay rate $\kappa/g$ increases from 0 to 0.15. The interference fringes vanish first (they involve high-order coherences), while the two Gaussian lobes persist — the state decoheres into a classical mixture. **Right:** Wigner negativity volume $\delta$ tracking the continuous loss of non-classicality.

### Dressed-State Avoided Crossing

<p align="center">
  <img src="../animations/anim_avoided_crossing.gif" width="750">
</p>

The JC dressed states $|n, \pm\rangle$ are the exact eigenstates of the coupled atom-cavity system. As the atom-cavity detuning $\Delta = \omega_a - \omega_c$ is swept, the bare-state energies (dashed grey) would cross, but the JC interaction opens an avoided crossing with a gap of $2g\sqrt{n+1}$. This animation sweeps $\Delta/g$ from $-10$ to $+10$ for manifolds $n = 0, 1, 5, 10$ simultaneously, making the $\sqrt{n+1}$ scaling of the vacuum Rabi splitting directly visible. At large detuning the dressed states approach the bare (uncoupled) states; on resonance the hybridization is maximal.

### Photon Number Distribution Dynamics

<p align="center">
  <img src="../animations/anim_photon_number.gif" width="750">
</p>

Time evolution of the intracavity photon number distribution $P(n, t) = \langle n|\rho_\text{field}(t)|n\rangle$ for an initial coherent state with $\bar{n} = 10$. At $t = 0$ the distribution is Poissonian (red dashed envelope). The two branches of the JC cat have the *same* amplitude $|\beta| = \sqrt{\bar{n}}$ and opposite phases ($\beta = \pm i\sqrt{\bar{n}}$), so their Poisson envelopes coincide; what $P(n,t)$ shows at $t = t_r/2$ is their **parity interference**: an odd-$n$ comb ($P(\text{odd}) = 0.82$, $\langle(-1)^{\hat n}\rangle = -0.64$), the photon-number signature of a (nearly odd) cat state. At the first revival the comb washes out and a broadened single-hump distribution returns; it never fully returns to Poissonian because of the anharmonic $\sqrt{n+1}$ Rabi spectrum.

### Bloch Sphere Trajectory

<p align="center">
  <img src="../animations/anim_bloch_sphere.gif" width="750">
</p>

The reduced atomic state $\rho\_{\text{atom}} = \mathrm{Tr}\_{\text{field}}[\rho]$ traces a trajectory inside the Bloch sphere. A pure atomic state sits on the surface ($|\mathbf{r}| = 1$); entanglement with the field pulls the Bloch vector toward the center ($|\mathbf{r}| \to 0$, maximally mixed). **Left:** 3D Bloch sphere trajectory color-coded by time. The atom starts at the excited state (red dot, north pole) and spirals inward during the collapse (entanglement with the field). Near $t = t_r/2$ (gold star) it returns close to the surface ($|\mathbf{r}| = 0.96$): the atom is nearly pure again because it has **disentangled** from the field, which is then in a cat state. Around the first revival (green triangle) the Bloch vector plunges back inside as the branches rephase and re-entangle. **Center:** Bloch vector length $|\mathbf{r}|(t)$. **Right:** Von Neumann entropy $S(\rho\_{\text{atom}})$.

### Wigner vs Husimi Q-Function

<p align="center">
  <img src="../animations/anim_q_vs_wigner.gif" width="750">
</p>

Side-by-side comparison of the Wigner function $W(x,p)$ (top row) and Husimi Q-function $Q(\alpha) = \langle\alpha|\rho|\alpha\rangle/\pi$ (bottom row) during JC evolution with $\bar{n} = 10$. The Wigner function takes negative values — the negativity volume $\delta$ (shown in each panel) quantifies non-classicality. The Husimi Q is a Gaussian-smoothed Wigner function ($Q = W * G_\text{vacuum}$) and is non-negative by construction ($Q \geq 0$ always). At $t = t_r/2$ the Wigner function resolves the interference fringes sharply, while the Q-function shows only two smooth lobes. This demonstrates that Wigner negativity, not Q-function structure, is the proper witness of quantum coherence in phase space.

### Cat-State Decoherence Sweep

<p align="center">
  <img src="../animations/anim_phase_diagram.gif" width="800">
</p>

The cat state at $t = t_r/2$ is progressively destroyed as cavity decay $\kappa/g$ increases from 0 to 0.12, with $\bar{n} = 10$ held fixed. **Left:** Wigner function $W(x,p)$ of the reduced cavity field. At $\kappa = 0$, the full interference pattern is visible between the two coherent lobes — deep negative fringes certifying a macroscopic quantum superposition. As $\kappa$ increases, the fringes wash out first (they are encoded in high-order off-diagonal elements $\langle n|\rho|n + 2k\rangle$ with $k \gg 1$), while the two classical lobes persist longer. By $\kappa/g \approx 0.04$ the negativity volume has dropped below our observability threshold $\delta = 0.05$ (and to $\approx 0.002$ by $\kappa/g = 0.1$) — the cat has decohered into an essentially classical mixture. **Right top:** Wigner negativity volume $\delta(\kappa/g)$ tracing the quantitative loss of non-classicality; the horizontal grey line marks $\delta = 0.05$, our operational threshold for observability. **Right bottom:** Field-state purity $\mathrm{Tr}[\rho_\text{field}^2]$, which drops from ~0.96 (near-pure cat state) to ~0.38 (highly mixed), confirming that decoherence (fringe erasure) proceeds much faster than energy dissipation (photon loss).

### Photon Blockade Transition

<p align="center">
  <img src="../animations/anim_g2_blockade.gif" width="800">
</p>

The transition from weak to strong coupling as $g/\kappa$ is swept from 0.1 to 12. **Left:** Energy-level diagram of the first three JC manifolds ($n = 0, 1, 2$) in units of $\hbar g$, with linewidth bands (shaded) that visibly shrink as $\kappa/g$ decreases. A red arrow marks the coherent drive; the red ✗ marks the blocked second-photon transition once the anharmonic splitting exceeds the linewidth. **Right:** Cavity transmission spectrum $\langle n \rangle(\Delta)$ computed from the steady-state Lindblad equation at each $g/\kappa$. At weak coupling ($g \ll \kappa$), the spectrum is a single Lorentzian centered at $\Delta = 0$ — the atom is too weakly coupled to modify the cavity response. As $g$ crosses $\kappa$, the peak broadens and flattens. At $g/\kappa \gtrsim 2$, the spectrum splits into two resolved peaks at $\Delta = \pm g$ — the **vacuum Rabi doublet**, the spectroscopic signature of strong coupling. The splitting grows as $2g$ (red annotation), directly mirroring the dressed-state gap in the energy ladder. Numerical readouts track $g^{(2)}(0)$ with the drive resonant with the lower polariton ($\Delta = \omega_c - \omega_L = g$), dropping from ~1 to deep antibunching (0.03 at $g/\kappa = 10$), and $\langle n \rangle$ at that detuning. Driving on *bare-cavity* resonance ($\Delta = 0$) in strong coupling gives the opposite: the laser is detuned from both polaritons, $\langle n\rangle \sim 10^{-8}$ and $g^{(2)}(0) \gg 1$ (bunching), so that is not the blockade signature.

---

---

## Static Figures

### Dressed-State Avoided Crossing and $\sqrt{n+1}$ Scaling

<p align="center">
  <img src="../figures/fig_avoided_crossing.png" width="800">
</p>

**Left:** Dressed-state energy eigenvalues $E_{n,\pm}$ as a function of detuning $\Delta/g$ for photon manifolds $n = 0$ (blue), $1$ (orange), $5$ (green), $10$ (red). At resonance ($\Delta = 0$), each manifold exhibits an avoided crossing with splitting $\Omega_n(0) = 2g\sqrt{n+1}$, indicated by colored arrows. **Right:** On-resonance splitting $\Omega_n(0)/g$ vs photon number $n$ (black dots), overlaid with the analytic curve $2\sqrt{n+1}$ (dashed pink). The $\sqrt{n+1}$ dependence is the quantum-mechanical fingerprint of the quantized field: a classical drive would produce a splitting independent of intensity.

### Photon Number Distribution at Six Key Times

<p align="center">
  <img src="../figures/fig_photon_number_evolution.png" width="800">
</p>

Snapshots of the intracavity photon number distribution $P(n,t)$ at six characteristic times during JC evolution ($\bar{n} = 10$, $\Delta = 0$). **(a)** $t = 0$: initial Poissonian distribution. **(b)** $t = 0.5\,t_c$: early Rabi oscillations, still approximately unimodal. **(c)** $t = 2\,t_c$: collapse onset, broadening as different Fock components oscillate at incommensurate $\sqrt{n+1}$ frequencies. **(d)** $t = t_r/2$ (cat state): odd-$n$ parity comb (the two cat branches have equal $|\beta|$, so the signature is interference between them, not a bimodal split). **(e)** $t = 0.75\,t_r$: partial recombination. **(f)** $t = t_r$ (first revival): partially recovered unimodal shape, broader than the initial Poissonian.

### Bloch Sphere Dynamics of the Reduced Atomic State

<p align="center">
  <img src="../figures/fig_bloch_sphere_trajectory.png" width="800">
</p>

**Left:** 3D trajectory of the reduced atomic Bloch vector $\mathbf{r} = (\mathrm{Tr}[\rho_\text{atom}\sigma_x],\, \mathrm{Tr}[\rho_\text{atom}\sigma_y],\, \mathrm{Tr}[\rho_\text{atom}\sigma_z])$ during one full collapse-revival cycle. The atom starts at the excited state (red dot, $|\mathbf{r}| = 1$), spirals inward during the collapse, returns to $|\mathbf{r}| = 0.96$ near the surface at $t = t_r/2$ (gold star: atom nearly pure, disentangled from the field cat), and plunges back inside around $t = t_r$ (green triangle: re-entangled). **Center:** Bloch vector length $|\mathbf{r}|(t)$. **Right:** Von Neumann entropy $S(\rho_\text{atom})$. Vertical dashed lines mark $t_c$ (blue), $t_r/2$ (pink), and $t_r$ (purple).

### Wigner vs Husimi Q-Function Comparison

<p align="center">
  <img src="../figures/fig_q_vs_wigner.png" width="800">
</p>

Side-by-side snapshots at six characteristic times. Each Wigner panel is annotated with the negativity volume $\delta$. At $t = t_r/2$, the Wigner function exhibits oscillatory fringes with $\delta = 0.85$, while the Q-function shows only two smooth, positive peaks. At $t = t_r$, residual negativity ($\delta = 0.31$) reflects imperfect refactorization.

### Cat-State Survival Phase Diagram

<p align="center">
  <img src="../figures/fig_phase_combined.png" width="900">
</p>

Systematic 2D parameter sweep (16 × 16 = 256 independent Lindblad simulations) mapping cat-state survival in the $(\bar{n},\, \kappa/g)$ plane at $t = t_r/2$. **(a)** Wigner negativity volume $\delta(\bar{n}, \kappa/g)$. The cyan contour marks $\delta = 0.05$ — the practical boundary below which cat-state interference fringes are unobservable. At $\kappa = 0$, negativity grows with $\bar{n}$ and saturates near $\delta \approx 0.88$ for $\bar{n} \gtrsim 14$. Any nonzero $\kappa$ degrades the cat: the $\delta = 0.05$ boundary follows $\kappa_c/g \approx 0.46\,\bar{n}^{-1.14}$ (fit over $\bar{n} = 4$–$25$), between the naive $1/\bar{n}$ and the $\bar{n}^{-3/2}$ estimate of Sec. 3 below. **(b)** Cat fidelity $F_\text{cat} = \max_{\beta,\theta}\langle\text{cat}|\rho_\text{field}|\text{cat}\rangle$ with $|\text{cat}\rangle \propto |\beta\rangle + e^{i\theta}|-\beta\rangle$ and $\beta$ near $i\sqrt{\bar{n}}$. The JC cat lives on the $p$ axis and is approximately *odd*, so comparing with the even cat $|\alpha\rangle + |-\alpha\rangle$ on the real axis would give $F \approx 0$ for every $\bar{n} \gtrsim 3$. Without dissipation $F_\text{cat} \approx 0.78$ for every $\bar{n}$ from 2 to 25: the JC interaction produces a cat of essentially fixed quality, and the $F = 0.5$ contour (cyan) tracks the $\delta = 0.05$ contour.

### Cat-State Survival: Parameter Slices

<p align="center">
  <img src="../figures/fig_phase_slices.png" width="900">
</p>

1D slices through the phase diagram. **Left:** $\delta$ vs $\bar{n}$ at fixed $\kappa/g$. Without dissipation ($\kappa = 0$, blue), negativity grows with $\bar{n}$ and saturates near $\delta \approx 0.88$ for $\bar{n} \gtrsim 14$. Modest decay ($\kappa/g = 0.03$, orange) limits the range with $\delta > 0.05$ to $\bar{n} \lesssim 10$. The horizontal grey line marks $\delta = 0.05$. **Right:** $\delta$ vs $\kappa/g$ at fixed $\bar{n}$. The decay is approximately exponential in $\kappa$, with a decay constant that shrinks roughly as $\bar{n}^{-1.1}$ — the enhanced decoherence rate of macroscopic superpositions.

### Photon Blockade: $g^{(2)}(0)$ vs Coupling Strength

<p align="center">
  <img src="../figures/fig_g2_vs_coupling.png" width="800">
</p>

Equal-time second-order coherence $g^{(2)}(0)$ of the intracavity field as a function of the vacuum Rabi coupling $g/\kappa$, with the drive resonant with the lower polariton ($\Delta = \omega_c - \omega_L = g$), for several drive amplitudes $\varepsilon/\kappa$. At weak coupling ($g \ll \kappa$) the cavity responds linearly and $g^{(2)}(0) \to 1$ (a coherently driven damped cavity is in a coherent state). As $g$ exceeds $\kappa$ (red dashed line) the anharmonic JC ladder blocks the second photon and $g^{(2)}(0)$ falls steadily (0.11 at $g/\kappa = 5$, 0.03 at $g/\kappa = 10$ for $\varepsilon/\kappa = 0.05$); weaker drives give deeper antibunching. The grey dashed line shows the same sweep with the drive on *bare-cavity* resonance ($\Delta = 0$): the laser is then detuned from both polaritons while the two-photon transition $|0\rangle \to |2,\pm\rangle$ is only $g/\sqrt{2}$ away, so $\langle n\rangle$ collapses to $\sim 10^{-8}$ and $g^{(2)}(0)$ climbs to $\sim 10^7$ — strong *bunching*, not blockade. **Bottom:** mean intracavity photon number $\langle n \rangle$ at $\Delta = g$.

### Photon Statistics vs Drive Strength

<p align="center">
  <img src="../figures/fig_g2_vs_drive.png" width="800">
</p>

$g^{(2)}(0)$ as a function of drive strength $\varepsilon/\kappa$ for several values of $g/\kappa$, drive on the polariton ($\Delta = g$). At weak drive the blockade sets $g^{(2)}(0) \approx 0.5$, 0.10 and 0.03 for $g/\kappa = 2$, 5 and 10. As $\varepsilon/\kappa$ grows past $\sim 0.3$ the drive overwhelms the anharmonicity and $g^{(2)}(0)$ rises monotonically toward the coherent-state value 1. Stronger coupling ($g/\kappa = 10$, green) maintains antibunching to higher drive powers.

### Photon Blockade Spectrum

<p align="center">
  <img src="../figures/fig_g2_blockade_spectrum.png" width="800">
</p>

**Top:** $g^{(2)}(0)$ as a function of laser-cavity detuning $\Delta/g$ for $g/\kappa = 5$, $\varepsilon/\kappa = 0.05$ (log scale). The antibunching dips ($g^{(2)}(0) \approx 0.1$) occur *at* the polariton resonances $\Delta = \pm g$ (red dashed lines): the first photon is absorbed resonantly, the second is blocked by the anharmonic ladder. Between them, the two-photon resonances at $\Delta = \pm g/\sqrt{2}$ (orange dashed lines), where $2\omega_L$ matches $E_{2,\pm}$, produce strong bunching ($g^{(2)}(0) \sim 10^2$), and on bare resonance $\Delta = 0$ the drive is far from every single-photon transition, $\langle n\rangle$ is tiny and $g^{(2)}(0)$ reaches $\sim 10^7$. **Bottom:** vacuum Rabi doublet in the cavity transmission spectrum $\langle n \rangle(\Delta)$ — the splitting of $2g$ is the spectroscopic signature of strong coupling.

### $g^{(2)}$ Combined Summary

<p align="center">
  <img src="../figures/fig_g2_combined.png" width="900">
</p>

Four-panel summary of photon blockade physics (drive on the polariton, $\Delta = g$, unless stated): **(a)** blockade transition vs $g/\kappa$, with the bare-resonance ($\Delta = 0$) bunching curve for contrast, **(b)** quantum-to-classical crossover vs drive strength, **(c)** blockade spectrum, **(d)** vacuum Rabi splitting in transmission.

### Cat-State Detail

<p align="center">
  <img src="../paper/figures/fig_cat_state_detail.png" width="750">
</p>

Cross-section $W(x, p{=}0)$ through the Wigner function at the cat-state time $t = t_r/2$ ($\delta = 0.85$, purity 0.96, $F_\text{cat} = 0.78$). Deep negative fringes reach $W \approx -0.25$, with fringe spacing $\pi / \sqrt{2\bar{n}} \approx 0.70$.

### Entanglement Across Field States

<p align="center">
  <img src="../paper/figures/fig_entanglement_comparison.png" width="800">
</p>

Reduced-atom entropy $S(\rho_\text{atom})$ (solid) and logarithmic negativity $E_\mathcal{N}$ (dashed) for four initial fields with $\bar{n} = 10$. **Coherent** (blue): collapse to ~1 bit, then a minimum $S \approx 0.14$ bit at $t_r/2$ (atom–field disentanglement; the field is the cat), then re-entanglement toward $t_r$. **Thermal** (red): $S$ saturates at 1 bit but $E_\mathcal{N} \approx 0.2$ — the global state is mixed and most of $S$ is classical correlation. **Squeezed vacuum** (green): partial dips only. **Fock $|10\rangle$** (panel b): strictly periodic, inversion period $\pi/(g\sqrt{11})$ and entropy period $\pi/(2g\sqrt{11})$. Thermal and squeezed states use $N_\text{cav} = 150$ ($N_\text{cav} = 50$ would truncate 1–3% of their weight and lower $\langle n\rangle$ to 9.6 and 9.2).

### Entropy Scaling with Photon Number

<p align="center">
  <img src="../paper/figures/fig_entropy_nbar_scaling.png" width="800">
</p>

The collapse time $t_c \sim 1/g$ is independent of $\bar{n}$, while the half-revival time $\pi\sqrt{\bar{n}}/g$ (dotted) and $t_r = 2\pi\sqrt{\bar{n}}/g$ (dashed) scale as $\sqrt{\bar{n}}$. The entropy minimum at $t_r/2$ (disentanglement) deepens with $\bar{n}$ and is cleanly resolved for $\bar{n} \gtrsim 9$.

### Dissipative Entanglement

<p align="center">
  <img src="../paper/figures/fig_dissipative_entanglement.png" width="800">
</p>

Inversion, reduced entropy and logarithmic negativity $E_\mathcal{N}$ under cavity loss. For $\kappa > 0$ the global state is mixed, so $S$ no longer measures entanglement: for $\kappa/g \geq 0.05$ it saturates near 1 bit at late times while $E_\mathcal{N}$ decays to zero. The inversion revival is gone by $\kappa/g = 0.05$, whereas $E_\mathcal{N}(t_r/2) \approx 0.35$ is almost unchanged up to $\kappa/g = 0.1$.

### Decoherence Table

| $\kappa/g$ | Wigner negativity $\delta$ | Field purity | Cat state visible? |
|-----------|--------------------------|--------------|-------------------|
| 0.00      | 0.851                    | 0.960        | Yes               |
| 0.02      | 0.143                    | 0.460        | Marginal          |
| 0.05      | 0.022                    | 0.418        | No                |
| 0.10      | 0.002                    | 0.386        | No                |

$\delta = \int|W|\,dx\,dp - 1$; values converged to three digits between $N_\text{cav} = 35$ and 50.

### Coherent vs Thermal Inversion

<p align="center">
  <img src="../paper/figures/jaynes_cummings_comparison.png" width="800">
</p>

Atomic inversion $\langle\sigma_z\rangle(t) = \sum_n P(n)\cos(2g\sqrt{n{+}1}\,t)$ for coherent (Poisson) and thermal (Bose-Einstein) initial fields at $\bar{n} = 4, 9, 14, 19, 24$; time in units of $1/g$.

### Mollow Triplet

<p align="center">
  <img src="../paper/figures/mollow_triplet_driving_strength.png" width="800">
</p>

Resonance fluorescence spectrum via the quantum regression theorem. Sidebands at $\pm\Omega$ with HWHM $= 3\gamma/4$, peak ratio 3:1.

### Static Wigner Functions for Fock States

<p align="center">
  <img src="../paper/figures/wigner_fock_combined.png" width="800">
</p>

$W_n(x,p) = \frac{(-1)^n}{\pi} L_n(2r^2) e^{-r^2}$ for Fock states $|n\rangle$. Ring structure and negativity grow with $n$.

---

