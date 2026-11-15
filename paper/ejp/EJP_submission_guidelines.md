# Submitting to the European Journal of Physics (IOP Publishing)

**Two candidate manuscripts now live here.** The newer one is the stronger submission:

- `jc_beyond_linear.tex` / `.pdf` (draft of 2026-11-13, proofread and corrected on
  2026-11-15): *Rotation and shear in the Jaynes-Cummings model*. About 4500 words of body
  text, an abstract of about 210 words, five figures, two tables, 34 references in citation
  order. It is built around a derivation: the second order of the expansion of the Rabi
  frequencies shears the two field components by an amount equal to their rotation angle, and
  closed forms follow for the revival heights, the parity at t_r/2, the cat fidelity (0.786),
  the squeezing, the Wigner negativity (0.869), the residual entanglement and the decay under
  loss, each tested against the exact solution. `jc_beyond_linear_supp.tex` / `.pdf` is its
  supplementary material (10 pages): derivations, numerical methods and their measured errors,
  limits, other initial states, notes on the exercises. Both are hand-written (not generated
  by `tools/`); figures come from `simulations/sim_rotation_shear.py`,
  `simulations/sim_numerics_checks.py` and `simulations/sim_phase_diagram.py`. The supplement
  takes the main text's equation, figure and table numbers from `jc_main_refs.tex`: after any
  renumbering run `python3 tools/sync_supp_refs.py`, then rebuild the supplement.
  Proofreading pass of 2026-11-15 (five independent checks: theory, recomputed numbers,
  references and attribution, language, figures; each finding was re-verified before it was
  applied). Its main result concerns priority:
  - **Gea-Banacloche, Phys. Rev. A 44, 5913 (1991), Appendix B already contains the
    second-order analysis** (read from the page images): the state with the quadratic phase
    (B4); the quadrature variances of each component (B11)-(B13), which are the eigenvalues of
    our V(s) (6.26 dB at t_r/2); the Q function (B16)-(B18), whose peak is (1+phi^2/4)^(-1/2)
    (0.786 at t_r/2); the value 0.55 for the overlap at the first revival; and the analogy with
    a Kerr medium. The abstract, introduction, sections 2-5 and the conclusions now say so.
  - Also credited now: the quadratic phase as the source of "phase spreading" (Auffeves et al
    2003, Eq. (2); Bild et al 2023); the loss exponent (Meunier et al 2006, Eq. (45)) and the
    rate 2 kappa nbar sin^2(phi) (Auffeves et al 2003); the coherence decay rate between two
    coherent states (Walls and Milburn 1985); the leading parity -cos(pi nbar) (Assemat et al
    2019); the explicit envelope (Pavlik 2024) and its widths (Feranchuk and Leonov 2009); the
    Gaussian loss of purity (Gea-Banacloche 1992, Phys. Rev. A 46, 7307); a numerical cat
    fidelity of 0.79 (Munhoz and Vidiella-Barranco 2005).
  - Still not found anywhere, and presented as new: reading the second order as a shear equal
    to the rotation angle; the parity amplitude 0.733 and phase 0.502 at t_r/2; the negativity
    limit 0.869; the full 1/nbar coefficient (4+pi^2)/32 of the impurity; the correction j(k);
    the remark that the cat fidelity does not improve with nbar (contrary to Larson and
    Mavrogordatos 2021).
  - Corrections of substance: the closed forms are values at the nominal times k t_r, whereas
    the kth revival peaks near 2 pi k sqrt(nbar+1)/g and is 3% (nbar = 10) or 1% (nbar = 40)
    higher; the bending of the components is of relative order nbar^(-1/2) and is present
    already for a purely quadratic phase; the envelope is accurate to 1/nbar only at the centre
    of a bump (nbar^(-1/2) on the flanks); the impurity and the negativity do not follow a
    1/nbar law below nbar of several hundred; table 2 uses cats |beta> + e^(i theta)|-beta>
    (the free two-coherent-state optimum is 0.802, 0.793, 0.788 at nbar = 10, 25, 100, above
    the limit); converged negativities 0.8505, 0.8805, 0.8672 at nbar = 10, 25, 100;
    kappa_c nbar^(3/2)/g = 1.15 +- 0.03 for 12 <= nbar <= 25; truncated squeezed vacuum
    <n> = 8.23; the loss sweep has 13 to 8 grid points per fringe, not 14.
  - Consequence for the venue: with the second order already in print the paper is a
    pedagogical synthesis with a few new closed forms. That suits Eur. J. Phys. or
    Am. J. Phys.; for Phys. Scr. or J. Phys. B it would need the detuned extension below, and
    the teaching section and exercises would go.
  - Still to do by hand: read Gea-Banacloche, Opt. Commun. 88, 531 (1992) and Haroche and
    Raimond Sec. 7.4 before claiming anything in the detuned case; confirm from the full texts
    the n^(3/4) squeezing time (Woods and Gea-Banacloche 1993) and the correspondence with
    Kennedy et al's Bloch vector; fill in the funding, conflict-of-interest and ORCID fields
    and the journal name in the class header; cut a new Zenodo release and update the DOI
    (items 1 and 3-7 of section 3 below apply unchanged).
  First literature check, 2026-11-13: superseded by the list above.
  Second literature check, later on 2026-11-13 (arXiv full texts read for the Paris-group and
  circuit-QED papers; changes to the list above):
  - The decay of the coherence under loss for ANY time, exponent
    (2 kappa n^(3/2)/g)(phi - sin(2 phi)/2), is Eq. (45) of Meunier, Le Diffon, Ruef,
    Degiovanni and Raimond, Phys. Rev. A 74, 033802 (2006) (arXiv:quant-ph/0603191; their phi
    is 2 phi and their g is 2 g here). Our visibility formula is its value at t_r/2, so it is a
    re-derivation: the draft now cites this paper next to Gea-Banacloche 1993. Check Eq. (45)
    yourself. The large-field, n-independent t^3 law is in Bonzom, Bouzidi and Degiovanni,
    Eur. Phys. J. D 47, 133 (2008) (arXiv:0704.3489).
  - The second-order phase exp[+-i phi (n - nbar)^2 / 4 nbar], i.e. shear = rotation angle, is
    written down as the source of "phase spreading" in Auffeves et al, Phys. Rev. Lett. 91,
    230405 (2003), Eq. (2), and in Haroche's College de France lecture of 10 Dec 2002,
    Eq. (8-17); neither extracts a shear parameter, a fidelity or a negativity from it. Cited
    where the expansion is introduced (done 2026-11-15).
  - The cat fidelity at t_r/2 is in print as a NUMBER: 0.7872 for alpha = 5 in Table 2 of
    Munhoz and Vidiella-Barranco, J. Mod. Opt. 52, 1557 (2005) (now cited). The closed form
    (1 + pi^2/16)^(-1/2) was not found.
  - Bild et al, Science 380, 274 (2023) make cats by resonant Jaynes-Cummings interaction in a
    mechanical resonator, describe the components as "distorted (and squeezed)" at the cat
    time, quote the revival bound (1 + pi^2)^(-1/4) from Eberly et al 1980, and report
    fidelities of 76% (to the JC state) and 66% (to an ideal cat): the best experimental
    motivation for a closed-form fidelity limit.
- APS version of the same article (2026-11-15): `paper/aps/jc_beyond_linear_aps.tex` / `.pdf`
  (REVTeX 4.2 `reprint`, two columns, 10 pages) and `jc_beyond_linear_aps_supp.tex` / `.pdf`
  (Supplemental Material, one column, 10 pages). Both are GENERATED from the IOP sources above by
  `tools/make_aps_revtex.py`; never edit them by hand. The physics, numbers and wording are those
  of the source. The script changes only the form: class and title block, `figure*` and
  `ruledtabular` floats, Physical Review cross-references (Eq., Fig., Sec., Ref., Roman table and
  section numbers), the supplement cited through the reference list, American spelling and
  serial commas, nine displayed equations reset on two lines for the 3.4 in column (numbering
  unchanged), the reference list in Physical Review format in citation order, and acknowledgments
  followed by a data-availability section. Every rewrite is asserted, so the script stops if the
  source changes where it rewrites (an equation it re-lays out, a phrase with a serial comma, a
  reference): fix the table it names in the script and run again.
  The same script also writes a research-journal cut, `jc_beyond_linear_aps_research.tex` /
  `.pdf` (10 pages) and `jc_beyond_linear_aps_research_supp.tex` / `.pdf` (9 pages). It removes
  the section "Using this material in teaching" with its seven exercises, the supplement's
  "Notes on the exercises", the closing paragraph on teaching, the keyword "computational physics
  education", the abstract's sentence on the intended readers, the sentence listing assumed
  background, and the pointers to exercises 2 and 7; three sentences that spoke of students are
  reworded neutrally (the list is `CUT_*` in the script). Equations, figures, tables, numbers
  and references are unchanged; Conclusions becomes Sec. IX. The physics still lacks the detuned
  extension below, which a Physical Review A submission would need.
- Exploratory extension (2026-11-13), `simulations/sim_detuned_shear.py` ->
  `paper/figures/fig_detuned_shear`; in no manuscript yet. With detuning Delta,
  n_crit = Delta^2/4g^2, sin(theta) = (1 + n_crit/n)^(-1/2): rotation phi = g^2 t/Omega,
  shear s = phi n/(n + n_crit) (resonant value phi, dispersive Kerr limit phi n/n_crit);
  detecting the atom in |g> heralds a balanced cat with probability sin^2(theta)/2 and
  fidelity (1 + s^2/4)^(-1/2); with loss the fidelity is that times (1 + exp(-Gamma))/2,
  Gamma = (2 kappa n^(3/2)/(g sin theta))(phi - sin(2 phi)/2). All checked against the exact
  amplitudes or the master equation (for Gamma < 1, agreement to about 1% at nbar >= 25 and
  2% at nbar = 10).
  Consequences for a cat of fixed size D = 2 sqrt(n) sin(phi): a larger resonant field removes
  the shear (F -> 1 - D^2/32n) and lowers Gamma from pi kappa D^3/8g (the half-revival cat)
  to kappa D^3/6g; detuning helps only when the field is capped and the loss weak
  (kappa D n/g < 3/2 for small phi), and always costs heralding probability. There is no
  universal optimal detuning.
  Prior work: the off-resonant revival amplitude is Feranchuk and Leonov 2009, Eqs. (20)-(21);
  the dispersive Kerr term g^4/Delta^3 is Boissonneault, Gambetta and Blais, Phys. Rev. A 79,
  013819 (2009); "resonance is optimal for large cats" and the size limit
  n < (Omega_0 T_c/2 pi)^(2/3) are in Raimond and Haroche, Seminaire Poincare 2, 25 (2005);
  heralded cats near resonance were made by Deleglise et al, Nature 455, 510 (2008) (adiabatic
  entry with Ramsey pulses, not the sudden switch assumed here). NOT found: the shear law
  n/(n + n_crit), the heralding probability and fidelity, the detuned loss exponent, the
  trade-off rules. Unread and able to pre-empt the detuned two-branch formulas:
  Gea-Banacloche, Opt. Commun. 88, 531 (1992) (full text), Haroche and Raimond's book
  Sec. 7.4. The two-level formulas do not carry over to transmons in the dispersive regime.
  Fixed on 2026-11-15: `cat_fidelity` in that script stopped short of the optimum whenever a
  starting coordinate was exactly zero (phi = pi/2), because SciPy's default Nelder-Mead simplex
  is then flat in that direction. It now starts from an explicit simplex, and the script checks
  it against a differential-evolution search on every run. Only the phi = pi/2 numbers changed:
  heralded-cat fidelity 0.802 (was 0.780) at nbar = 10 on resonance, 0.936 (was 0.910) at
  nbar = 10 with Delta = 2 g sqrt(nbar), 0.788 (was 0.786) and 0.932 (was 0.929) at nbar = 100;
  the nbar = 10 points of panel (e) rose by up to 0.02. At finite nbar the best pair of coherent
  states sits slightly off the centres of the components, so the exact fidelity lies above the
  closed form and approaches it as 0.16/nbar. `fig_detuned_shear`,
  `paper/notes/detuned_derivation` and `paper/summary/jc_summary` were updated.
- `paper/qsol_cqed.tex` / `.pdf`, the REVTeX 4.2 (APS style) paper from which the three-lessons
  files below are generated, was proofread on 2026-11-15 by five independent checks (theory,
  recomputed numbers, references and experimental facts, language, figures) and corrected.
  Format: formal author name, no URL in the abstract, a Data availability statement, figures in
  the text in citation order (the 3D rendering of the static Wigner functions was dropped as a
  duplicate of the maps; 11 figures remain), REVTeX ruled tables, 46 references in citation
  order. Corrections of substance: the collapse-envelope carrier is Omega_0 sqrt(nbar+1) and
  the envelope exp[-g^2 t^2 nbar/2(nbar+1)]; the revival envelope peaks at 2 pi sqrt(nbar+1)/g;
  the mixing angle is defined on [0, pi]; the truncated squeezed vacuum at N_cav = 50 has
  <n> = 8.23 (9.15 only with a truncated squeeze operator); S(t_r/2) = 0.26 and 0.48 bit for
  kappa/g = 0.05 and 0.10; in the Paris experiments Omega_0/2pi = 50 kHz means g/2pi = 25 kHz
  here, with cavity damping times of 0.85 ms (2003) and 0.13 s (2007); the 1992 optical
  experiment used an atomic beam with one atom on average, not a single atom.
  The three-lessons files below were NOT regenerated after this pass. `tools/reframe_ejp.py`
  and `tools/make_ejp_paper.py` replace spans keyed on sentences that no longer exist in the
  source, and their own replacement texts repeat several statements corrected here (Paris
  parameters, "over two decades", "single cesium atom", "maximally entangled", collapse
  "complete by 2/g"). `tools/make_ejp_iopjournal.py` already reads the new source layout.
  Re-key and correct those two scripts before regenerating, or retire the three-lessons
  versions; until then their .tex and .pdf are as built on 2026-11-13.
- The earlier "three lessons" article, in three files:


- `qsol_cqed_ejp_paper.tex` / `.pdf`: the **submission candidate**, a "Paper" of 3999 words of
  body text, eight figures, two tables, in IOP's `iopjournal` class.
- `qsol_cqed_ejp.tex` / `.pdf`: the long version (about 9000 words) with the full theory
  section; useful as an arXiv companion or as source for supplementary material.
- `qsol_cqed_ejp_iopart_legacy.tex`: an old `iopart` copy with no generator. It was not
  updated in the November 2026 pass and is out of date; delete it or regenerate it by hand.

The journal facts below were checked on 2026-11-12 against the journal's "About" and author
guideline pages. Policies change: re-read https://publishingsupport.iopscience.iop.org/journals/european-journal-of-physics/
and https://iopscience.iop.org/journal/0143-0807 before submitting.

## 1. How the files are built

Do not edit the `.tex` files in this folder. They are generated:

```
paper/qsol_cqed.tex  --make_ejp_iopjournal.py-->  qsol_cqed_ejp.tex   (markup only)
                     --reframe_ejp.py-->          qsol_cqed_ejp.tex   (teaching framing, Sec. 7)
                     --make_ejp_paper.py-->       qsol_cqed_ejp_paper.tex  (condensed to ~4000 words)
```

Run the three scripts in `tools/` in that order from any directory, then
`tectonic qsol_cqed_ejp_paper.tex`. Each script checks that every `\ref` and `\cite` resolves
and `make_ejp_paper.py` prints the word count per section. Sentences wrapped in
`\revtexonly{...}` in the REVTeX source (they cite the four supplementary figures) are dropped.

- `iopjournal.cls` is not in TeX Live; a copy is in this folder and Overleaf's "IOP
  Publishing" template ships it together with the `orcid` icon needed by `\orcid{}`.
- Figures are included as vector PDF at `\textwidth` from `../figures/`. They are drawn at
  153 mm with 8 to 10 pt text, panels are lettered (a), (b), ..., and every series differs in
  line style or lightness as well as hue.
- No BibTeX: the bibliography is inline in IOP numeric style. Do not add a
  `\section*{References}` heading; the class prints one.

## 2. What the journal asks for

- **Scope.** EJP exists to maintain and improve the standard of taught physics in higher
  education and says it is not a vehicle for publishing fundamental research. A submission
  must say where its physics sits in the curriculum, how it advances teaching, why it is
  useful, and at what level (undergraduate or graduate); submissions that do not are not
  considered. "Successful original student projects (experimental, theoretical or
  computational)" are explicitly in scope, as are quantum optics and quantum mechanics.
- **Article type and length.** "Paper", typically no more than 4000 words; abstract up to 300
  words. Reviews (8000 to 10 000 words) are a different type, usually aimed at introducing a
  field.
- **Review.** Single-anonymous. Keep the author block; the `[anonymous]` class option is not
  needed.
- **First submission.** One PDF with the figures embedded, text at least 12 pt. Source files
  are requested at revision.
- **Figures.** Vector EPS/PDF preferred, text 8 to 12 pt at final size, width about 8.5 cm or
  15 cm, information not carried by colour alone.
- **Other.** ORCID required for the submitting author. Data availability statement required
  (code counts as data). Subscription-route publication is free; open access costs an APC
  unless your institution has an agreement with IOP.
- **Odds.** The journal's own figures for the last calendar year: about half of submissions
  are rejected without review, about a quarter are accepted, and the median time to a first
  decision after review is about eight weeks.

## 3. Where the manuscript stands against that

Done:

- Teaching framing in abstract and introduction: audience, prerequisites, three lessons.
- "Using this material in teaching" section: where it fits, four misconceptions, six
  exercises, an assessment task.
- Length 3999 words of body text (at the limit: anything added must be offset); abstract 279 words.
- Figures redrawn for print (Sec. 1); numbers in text, tables and figures agree with the
  scripts.
- Data statement and software citation in place.

Still to do before submitting:

1. **Cut a new release.** The code changed after Zenodo v1.1.0, so the data statements and
   reference [nguyen2026code] point at a version that does not produce these figures.
   `CITATION.cff` was bumped to version 1.2.0 on 2026-11-15. Push, create the GitHub release
   `v1.2.0`, wait for Zenodo to mint the version DOI, then run
   `python3 tools/set_release_doi.py 1.2.0 10.5281/zenodo.NNNNNNNN`: it writes the version and
   DOI into `CITATION.cff`, `README.md`, both manuscripts and the two generators, and prints
   the commands that rebuild the PDFs. (A release cannot contain its own version DOI, hence
   the two steps. The three-lessons EJP files keep v1.1.0, the version they were built from.)
2. **Cite the prior treatments** a referee will know, and say what this paper adds to each.
   The closest are Larson and Mavrogordatos, *The Jaynes-Cummings Model and Its Descendants*
   (IOP, 2021), Sec. 2.2; Knight, Gerry, Birrittella and Alsing, JOSA B 41, C91 (2024);
   Gerry and Knight, Am. J. Phys. 65, 964 (1997); Buzek, Moya-Cessa, Knight and Phoenix,
   Phys. Rev. A 45, 8190 (1992); Eiselt and Risken, Phys. Rev. A 43, 346 (1991); Scheel,
   Eisert, Knight and Plenio, J. Mod. Opt. 50, 881 (2003). Check each against the publisher
   page before citing.
3. **Review copy at 12 pt**: `\documentclass[12pt]{iopjournal}` if the class accepts it,
   otherwise submit the PDF as is and say so in the cover letter (the template states that
   formatting to house style is not required).
4. Fill in `\orcid{}`, confirm `\email{}`, and uncomment `\funding{}` and `\roles{}`.
5. Add a conflict-of-interest sentence if the submission system asks for one in the text.
6. Supplementary material: solutions to the six exercises, and one notebook that reproduces
   Figs. 2, 5 and 8. Upload figure files separately as `figure1.pdf` ... `figure8.pdf`
   (copy them from `../figures/` in the order they appear in the PDF).
7. Cover letter: state the level, the course it came from, the teaching gap, that every
   figure is reproducible from the archived code, and that the work is a student project.
