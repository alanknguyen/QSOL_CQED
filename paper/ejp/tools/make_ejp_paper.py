"""Derive the ~4000-word EJP 'Paper' version from the reframed EJP file.
Cuts: theory condensed to what the results use; criteria section folded into
two sentences; experiment section reduced to the parameter mapping; intro and
conclusions trimmed.  Abstract, methods, results, teaching section unchanged."""
import re, pathlib

SRC = pathlib.Path('/Users/alanknguyen/Documents/qsol_cqed/paper/ejp/qsol_cqed_ejp.tex')
OUT = SRC.with_name('qsol_cqed_ejp_paper.tex')
s = SRC.read_text()

def replace_span(start, end, new):
    global s
    assert s.count(start) == 1, start[:60]
    i = s.index(start); j = s.index(end, i); s = s[:i] + new + s[j:]

def rep(old, new):
    global s
    assert s.count(old) == 1, old[:60]
    s = s.replace(old, new)

# ---------------------------------------------------------------- header comment
rep("% European Journal of Physics (IOP Publishing) version of qsol_cqed.tex,",
    "% European Journal of Physics 'Paper' (~4000 words) version of qsol_cqed.tex,")

# ---------------------------------------------------------------- introduction
INTRO1 = r"""The Jaynes--Cummings (JC) model~\cite{jaynes1963}, a single two-level atom coupled to a single quantized cavity mode, is the simplest fully quantum description of light--matter interaction and the standard first encounter with it in a quantum-optics course. It predicts phenomena with no classical analog: vacuum Rabi splitting, the collapse and revival of Rabi oscillations~\cite{eberly1980}, and the dynamical generation of Schr\"{o}dinger cat states in the cavity field~\cite{brune1996}, all confirmed in microwave cavity QED~\cite{rempe1987, haroche2006} and in circuit QED with superconducting qubits~\cite{blais2021}.

"""
replace_span('Coherent manipulation of quantum two-level systems', 'While the JC model appears in every', INTRO1)
rep(r"""Section~\ref{sec:theory} supplies a self-contained derivation of the rest, and readers who know the JC model can proceed directly to Section~\ref{sec:methods}.""",
    r"""Section~\ref{sec:theory} collects the results we need and points to full derivations.""")
replace_span('The paper is organized as follows.', '\n\\section{Theoretical Framework}',
    r"""Section~\ref{sec:theory} summarizes the theory, Section~\ref{sec:methods} the computational methods, and Section~\ref{sec:results} the results; Section~\ref{sec:experiment} maps them onto experimental platforms, Section~\ref{sec:teaching} suggests classroom use, and Section~\ref{sec:conclusions} concludes.
""")

# ---------------------------------------------------------------- theory (condensed)
THEORY = r"""\section{Theoretical Framework}
\label{sec:theory}

\subsection{The Jaynes--Cummings model}
\label{sec:jcm}
A single cavity mode is a harmonic oscillator with ladder operators $a$, $a^\dagger$, $[a, a^\dagger] = 1$, and Fock states $|n\rangle$ with $a^\dagger|n\rangle = \sqrt{n+1}\,|n{+}1\rangle$. In the dipole and rotating-wave approximations, a two-level atom with transition frequency $\omega_a$ coupled to the mode of frequency $\omega_c$ is described by the Jaynes--Cummings Hamiltonian~\cite{jaynes1963, gerry2004}
\begin{equation}
H_{\text{JC}} = \hbar\omega_c\, a^\dagger a + \tfrac{1}{2}\hbar\omega_a\sigma_z + \hbar g\left(a^\dagger \sigma^- + a\,\sigma^+\right),
\label{eq:HJC}
\end{equation}
where $\sigma^\pm$ raise and lower the atom, $\sigma_z = |e\rangle\langle e| - |g\rangle\langle g|$, and $g$ is the vacuum Rabi coupling. $H_{\text{JC}}$ conserves the excitation number $a^\dagger a + \sigma^+\sigma^-$, so it decomposes into $2\times 2$ blocks spanned by $\{|e,n\rangle, |g,n{+}1\rangle\}$. With detuning $\Delta = \omega_a - \omega_c$ and $\Omega_0 = 2g$, each block has eigenvalues $E_\pm(n) = \hbar\omega_c(n + \tfrac{1}{2}) \pm \tfrac{1}{2}\hbar\,\Omega_n(\Delta)$ with the $n$-photon Rabi frequency $\Omega_n(\Delta) = \sqrt{\Delta^2 + \Omega_0^2(n+1)}$. Its eigenstates, the dressed states, are equal superpositions of $|e,n\rangle$ and $|g,n{+}1\rangle$ on resonance and approach the bare states in the dispersive limit $|\Delta| \gg \Omega_0\sqrt{n+1}$, where the coupling reduces to ac Stark shifts~\cite{haroche2006}.\label{sec:limits}

\subsection{Collapse and revival}
\label{sec:collapse_revival}
For the atom initially excited and the field in a Fock state $|n\rangle$, on resonance,
\begin{equation}
|\Psi(t)\rangle = \cos\!\left(\tfrac{1}{2}\Omega_n t\right)|e,n\rangle - i\sin\!\left(\tfrac{1}{2}\Omega_n t\right)|g,n{+}1\rangle, \qquad \Omega_n = \Omega_0\sqrt{n+1},
\label{eq:rabi_evolution}
\end{equation}
so the excited-state probability oscillates as $P_e(t) = \cos^2(\tfrac{1}{2}\Omega_n t)$; $n = 0$ is the vacuum Rabi oscillation. For a coherent field $|\alpha\rangle$ with Poissonian statistics $p(n) = e^{-\bar{n}}\bar{n}^n/n!$ and $\bar{n} = |\alpha|^2$, each Fock component evolves independently and
\begin{equation}
P_e(t) = \frac{1}{2} + \frac{1}{2}\sum_{n=0}^{\infty} p(n)\cos\!\left(\Omega_0\sqrt{n+1}\,t\right).
\label{eq:pe_cos_sum}
\end{equation}
The Rabi frequencies of the components spread by $\delta\Omega \approx \Omega_0/2$ across the Poisson width $\sqrt{\bar{n}}$, so the oscillations dephase on the collapse time $t_c \sim 1/g$, independent of $\bar{n}$; linearizing $\sqrt{n+1}$ about $\bar{n}$ gives the envelope $\tfrac{1}{2}\cos(\Omega_0\sqrt{\bar{n}}\,t)\,e^{-g^2t^2/2}$~\cite{eberly1980}. Throughout we use the operational definition $t_c \equiv 1/g$. Because the spectrum is discrete, adjacent components rephase when $(\Omega_0\sqrt{n+2} - \Omega_0\sqrt{n+1})\,t_r = 2\pi$, giving the revival time
\begin{equation}
t_r = \frac{2\pi\sqrt{\bar{n}}}{g},
\label{eq:t_revival}
\end{equation}
a direct signature of field quantization~\cite{eberly1980, haroche2006}. A thermal field, with Bose--Einstein statistics $p(n) = \bar{n}^n/(1+\bar{n})^{n+1}$ and width $\sim\bar{n}$, dephases within a Rabi period and shows no resolved revival; a Fock field never dephases. During the collapse the atom and field become strongly entangled. At half the revival time, $t_r/2 = \pi\sqrt{\bar{n}}/g$, the two branches of the atomic state reconverge, the atom and field approximately \emph{disentangle}, and the field is left in a Schr\"{o}dinger cat state, a superposition of two coherent states with opposite phases~\cite{gea-banacloche1990, gea-banacloche1991, phoenix1991}.

\subsection{Phase space, dissipation, and entanglement measures}
\label{sec:lindblad}
The Wigner function of the reduced field state $\rho_{\text{field}}$~\cite{wigner1932},
\begin{equation}
W(x, p) = \frac{1}{\pi\hbar}\int_{-\infty}^{\infty} \langle x - y | \rho_{\text{field}} | x + y \rangle\, e^{2ipy/\hbar}\, dy,
\label{eq:wigner_def}
\end{equation}
is the Gaussian $\pi^{-1}\exp[-(x - \sqrt{2}\,\text{Re}\,\alpha)^2 - (p - \sqrt{2}\,\text{Im}\,\alpha)^2]$ for a coherent state and takes negative values for Fock and cat states. Negativity is a sufficient but not a necessary witness of non-classicality~\cite{hudson1974}: the squeezed vacuum has a non-negative Gaussian $W$ yet a singular Glauber--Sudarshan $P$ function~\cite{glauber1963, sudarshan1963}. We quantify negativity by the volume~\cite{kenfack2004}
\begin{equation}
\delta = \int \big|W(x,p)\big|\, dx\, dp - 1,
\label{eq:negativity}
\end{equation}
which equals twice the integrated negative part of $W$ because $\int W\,dx\,dp = 1$.

Cavity photon loss at rate $\kappa$ and atomic decay at rate $\gamma$ enter through the Lindblad master equation~\cite{lindblad1976, wiseman2009}
\begin{equation}
\frac{d\rho}{dt} = -\frac{i}{\hbar}[H_{\text{JC}}, \rho] + \kappa\,\mathcal{D}[a]\rho + \gamma\,\mathcal{D}[\sigma^-]\rho, \qquad \mathcal{D}[L]\rho = L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L, \rho\}.
\label{eq:lindblad}
\end{equation}
The interference fringes of a superposition of two coherent states separated by $\Delta\alpha$ decay at the rate $\kappa|\Delta\alpha|^2/2$~\cite{zurek2003, haroche2006}; for the JC cat, whose components sit at $\pm i\sqrt{\bar{n}}$ (Sec.~\ref{sec:wigner_results}), this is $2\kappa\bar{n}$, far faster than the energy decay rate $\kappa$.

For a pure atom--field state the entanglement is measured by the von Neumann entropy of the reduced atom~\cite{nielsen2000},
\begin{equation}
S(\rho_{\text{atom}}) = -\text{Tr}\big[\rho_{\text{atom}} \log_2 \rho_{\text{atom}}\big],
\label{eq:entropy}
\end{equation}
which ranges from 0 to 1~bit and is complemented by the field purity $\mathcal{P} = \text{Tr}[\rho_{\text{field}}^2]$. If the global state is mixed, as for a thermal field or under dissipation, $S(\rho_{\text{atom}})$ counts classical correlations as well and can reach 1~bit for a weakly entangled state; we then use the logarithmic negativity~\cite{vidal2002, plenio2005}
\begin{equation}
E_{\mathcal{N}} = \log_2 \big\| \rho^{T_{\text{atom}}} \big\|_1,
\label{eq:logneg}
\end{equation}
an entanglement monotone that vanishes for separable states and equals 1 for a maximally entangled pair. Phoenix and Knight~\cite{phoenix1988, phoenix1991} first traced $S$ through the collapse--revival cycle; Sec.~\ref{sec:results} verifies that picture computationally.


"""
replace_span('\\section{Theoretical Framework}', '\\section{Computational Methods}', THEORY)

# ---------------------------------------------------------------- criteria section: fold into results, remove
rep(r"""The revival is not exact: higher-order terms in the photon-number distribution prevent perfect recurrence.""",
    r"""The revival is not exact: higher-order terms in the photon-number distribution prevent perfect recurrence. Negative $W$ is the phase-space face of the non-classicality that photon counting detects as sub-Poissonian statistics or antibunching~\cite{kimble1977}; the Glauber--Sudarshan $P$ function remains the definitive criterion, and the everywhere-positive Husimi $Q$ function is correspondingly less sensitive.""")
replace_span('\\section{Classical vs.\\ Quantum Light', '\\section{Connection to Experiment}', '')
s = s.replace('(Sec.~\\ref{sec:criteria})', '(Sec.~\\ref{sec:lindblad})')
assert 'sec:criteria' not in s

# ---------------------------------------------------------------- experiment (condensed)
EXPT = r"""\section{Connection to Experiment}
\label{sec:experiment}

Collapse and revival were first observed with Rydberg atoms crossing a superconducting microwave cavity~\cite{rempe1987}; the same platform later generated and detected cat states~\cite{brune1996, haroche2006}, trapped ions reproduced the full sequence with a vibrational mode playing the role of the field~\cite{meekhof1996}, and circuit QED now reaches $g/\kappa \sim 100$--$300$ with superconducting qubits~\cite{blais2021}, where cat states serve as error-corrected qubits~\cite{ofek2016}. Our dimensionless results map onto any platform by inserting the physical value of $g$. For microwave cavity QED ($g/2\pi \sim 50$~kHz, $\kappa/2\pi \sim 1$~Hz) the revival time for $\bar{n} = 10$ is $t_r = 2\pi\sqrt{10}/g \approx 60~\mu$s against a cavity lifetime $1/\kappa \approx 0.16$~s; for circuit QED ($g/2\pi \sim 100$--$300$~MHz, $\kappa/2\pi \sim 1$~MHz) it is $10$--$30$~ns against $160$~ns. The cat fringes, which decay at $2\kappa\bar{n}$, survive only where $\kappa/g \lesssim 0.01$--$0.03$, the threshold of Sec.~\ref{sec:decoherence}; this is met in both regimes but not in typical optical cavities.


"""
replace_span('\\section{Connection to Experiment}', '\\section{Using this material in teaching}', EXPT)

# ---------------------------------------------------------------- conclusions (trimmed)
CONCL = r"""\section{Conclusions}
\label{sec:conclusions}

We have presented a computational treatment of the Jaynes--Cummings model intended for teaching, in which the collapse and revival of Rabi oscillations, the phase-space dynamics of the cavity field, and the atom--field entanglement are followed together from reproducible code. Three lessons emerge. First, the Schr\"{o}dinger cat state that forms at half the revival time, with Wigner negativity $\delta = 0.85$, purity 0.96 and cat fidelity 0.78, coincides with a \emph{minimum} of the atom--field entanglement ($S \approx 0.14$~bit): it is a nearly pure state of the field alone, whereas during the preceding collapse the atom and field are maximally entangled and the reduced field is a two-branch mixture. Second, this half-revival disentanglement~\cite{gea-banacloche1990, phoenix1991} is unique to the coherent field, and for a thermal field a reduced entropy of one bit conceals a logarithmic negativity of only about 0.2: once the global state is mixed, the reduced entropy is not an entanglement measure. Third, cavity decay as weak as $\kappa/g = 0.02$ removes over 80\% of the Wigner negativity while leaving the residual atom--field negativity at $t_r/2$ nearly unchanged, separating the loss of the field's internal coherence from the loss of entanglement.

The physics is classical work; what the computational route adds is that the three observables can be placed side by side at negligible cost, so that each misconception listed in Sec.~\ref{sec:teaching} becomes a figure a student can regenerate and then break by changing a parameter.


"""
replace_span('\\section{Conclusions}', '\\ack{', CONCL)

# ---------------------------------------------------------------- prune uncited references
cites = set(k.strip() for m in re.findall(r'\\cite\{([^}]*)\}', s) for k in m.split(','))
dropped = []
for key in re.findall(r'\\bibitem\{([^}]*)\}', s):
    if key not in cites:
        s = re.sub(r'\\bibitem\{' + re.escape(key) + r'\}[^\n]*\n', '', s); dropped.append(key)

# ---------------------------------------------------------------- checks
labels = set(re.findall(r'\\label\{([^}]*)\}', s)); refs = set(re.findall(r'\\ref\{([^}]*)\}', s))
assert refs <= labels, refs - labels
bibkeys = set(re.findall(r'\\bibitem\{([^}]*)\}', s)); assert cites <= bibkeys, cites - bibkeys
envs = re.findall(r'\\begin\{([a-zA-Z*]+)\}', s); ends = re.findall(r'\\end\{([a-zA-Z*]+)\}', s)
assert all(envs.count(k) == ends.count(k) for k in set(envs) | set(ends))
t = re.sub(r'(?<!\\)%.*', '', s).replace('\\{', '').replace('\\}', '')
assert t.count('{') == t.count('}'), (t.count('{'), t.count('}'))
OUT.write_text(s)

# word counts per section (prose only)
body = s.split('\\end{abstract}')[1].split('\\begin{thebibliography}')[0]
body = re.sub(r'\\begin\{(figure|table)\}.*?\\end\{\1\}', ' ', body, flags=re.S)
def wc(t):
    t = re.sub(r'\$[^$]*\$', 'X', t); t = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?', ' ', t)
    return len(re.sub(r'[{}~&\\]', ' ', t).split())
parts = re.split(r'(\\section\{[^}]*\})', body); total = 0
for i in range(1, len(parts), 2):
    n = wc(parts[i+1]); total += n; print(f"{n:6d}  {re.search(r'section\{([^}]*)\}', parts[i]).group(1)}")
print(f"{total:6d}  TOTAL prose words (excl. abstract, captions, tables, references)")
print(f"wrote {OUT.name}; {len(bibkeys)} references kept, dropped uncited: {dropped}")
