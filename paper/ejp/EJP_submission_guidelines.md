# Submitting to the European Journal of Physics (IOP Publishing)

`qsol_cqed_ejp.tex` is the corrected article reformatted for IOP Publishing's current
journal template, `iopjournal.cls` (2024; the class Overleaf's "IOP Publishing" template
ships). The physics content is identical to `../qsol_cqed.tex`; only the markup changed.
`qsol_cqed_ejp_iopart_legacy.tex` is the same article for the old 1996 `iopart` class,
kept only in case a journal still asks for it.
This note lists what still has to be done by hand, and what the journal will look for.
Always check the current author guidelines at https://publishingsupport.iopscience.iop.org/
and the journal page https://iopscience.iop.org/journal/0143-0807 before submitting; policies change.

## 1. Building the file

- `iopjournal.cls` is not in TeX Live but needs only standard packages (article, fancyhdr,
  xcolor, graphicx, hyperref), so it compiles anywhere the class file sits next to the
  `.tex`. A copy is in this folder; Overleaf's IOP template provides it too, together with
  the `orcid` icon needed by `\orcid{}`.
- Overleaf: create the project from the IOP Publishing template, replace its `main.tex`
  with `qsol_cqed_ejp.tex`, and upload the eight figures into a folder named `figures/`
  (the file has `\graphicspath{{figures/}{../figures/}}`). Compile with pdfLaTeX.
- The file was compiled locally with the real class: no errors, all references resolved.
- No BibTeX is needed: the bibliography is inline in IOP numeric style. Do **not** add a
  `\section*{References}` line above it; in this article-based class `thebibliography`
  prints the heading itself and you would get it twice.
- The EJP file contains Figures 1–8 only, each placed after the paragraph that first cites
  it. Figures 9–12 of the REVTeX version (semi-analytic inversion, the two static Wigner
  galleries, the Mollow triplet) were never cited in the text and are omitted here; they
  stay in `../qsol_cqed.tex` and in the repository, and can be offered as supplementary
  material via `\suppdata{}`.
- Fill in or confirm: the e-mail in `\email{}`, the two `\affil{}` lines, an ORCID via
  `\orcid{}` after the author name (requires the template's icon file), and optionally
  `\funding{}` and `\roles{}` (both left as comments).
- The running head (`\fancyhead`) carries journal placeholders `vv (yyyy) aaaaaa`; the
  journal fills those in. `\articletype{Paper}` prints the template's "Journal Name" banner
  and a Crossmark/Received/Revised margin mock-up; that is how the template looks and is
  not something to edit.
- Double-anonymous review: `\documentclass[anonymous]{iopjournal}` removes the author
  list, affiliations, e-mail, acknowledgments, funding, roles and data statement
  automatically. Check whether EJP requires or offers anonymous review and use the option
  for the reviewed copy if so; also remove the repository URL from the abstract in that
  case, since it identifies the author.

## 2. What EJP actually publishes, and how this article must be reframed

EJP is a physics-education journal. Referees ask "what does this teach, to whom, and why is
it better than what exists?", not "what is new physics?". The current text still carries the
research-paper framing ("goes beyond standard treatments", "original computational results").
Before submission:

- **Rewrite the abstract and introduction for a teaching audience.** State the level
  (advanced undergraduate / first-year graduate quantum optics), the prerequisites
  (second quantization, density matrices, a first exposure to the Jaynes–Cummings model),
  and the learning goals: (i) collapse and revival as dephasing and rephasing of a discrete
  spectrum; (ii) the difference between a two-branch *mixture* and a coherent *superposition*,
  read off from purity and Wigner negativity; (iii) why the reduced entropy is an
  entanglement measure only for globally pure states, using the thermal case and the
  logarithmic negativity as the counterexample; (iv) half-revival disentanglement and cat
  formation as one event seen through three observables.
- **Add a short section "Using this material in teaching"** (2–3 paragraphs): which figures
  map to which lecture, suggested exercises students can run with the repository
  (change $\bar n$, change $\kappa/g$, reproduce Table 2, verify the entropy period of the
  Fock case), and typical run times. EJP values exercises with solutions in supplementary
  material.
- **Present the code as the pedagogical instrument.** Point to the repository early, give the
  QuTiP version, and consider archiving a tagged release on Zenodo so the Data availability
  statement can cite a DOI rather than a moving GitHub URL.
- **Cut to the story.** Figures 1–8 carry the argument and are the only ones in the EJP
  file. The omitted Figs. 9–12 can go to supplementary material or the repository notebooks
  if a teaching use is written for them. A 6–8 figure paper of about 5000–6000 words is
  typical for EJP; the present body is about 7600 words, so the text still needs trimming,
  mainly in Sections 2 and 6.
- **Tone down claims of novelty.** Cite Gea-Banacloche (1990, 1991), Phoenix and Knight
  (1988, 1991), Bužek et al (1992) and Eiselt and Risken (1991) as the sources of the
  physics, and present the article's contribution as a verified, reproducible computational
  route through that physics with modern metrics (negativity, log-negativity, cat fidelity).
- **Keep the corrections visible.** The three-regime picture (entangled mixture during the
  collapse, disentangled pure cat at $t_r/2$, re-entanglement at $t_r$) is the pedagogical
  core; the common misreading that the cat coincides with maximal entanglement is exactly
  the kind of misconception EJP likes to see resolved with a figure.

## 3. IOP formatting rules the file already follows

- Sentence-case title; author, affiliations, e-mail and keywords via the template's
  `\author`, `\affil`, `\email` and `\keywords` commands, then the abstract (the template
  says at most about 300 words; the current abstract is about 265, so trim it slightly when reframing for EJP).
- Sections numbered by the class; back matter via `\ack{}` and `\data{}` (a Data
  availability statement is required by IOP journals).
- Figures: `\begin{figure}\centering\includegraphics...\caption{...}\label{...}\end{figure}`.
  Colour is free online; make sure every figure is legible in greyscale (Figs. 4 and 8 use
  colour and line style, which is fine). IOP asks that colour is not the only carrier of
  information.
- Tables: plain `tabular` with `\hline` rules, caption above, as in the template.
- References: numeric, in order of citation, IOP style "Author A B, Author C D and Author E F
  year *Journal* **volume** page"; article titles are not required. All 37 entries were
  converted; entry `walls2008` is not cited in the text and can be dropped.
- The template itself says formatting to the published style is not required; clarity is
  what matters. Do not spend time on layout beyond what the class gives you.

## 4. Things to prepare besides the manuscript

- Figure files: PDF for line art (available for all figures) or PNG at ≥300 dpi; IOP accepts
  EPS/PDF/PNG/TIFF. Upload figures separately at submission even if embedded.
- Supplementary material: the repository link, a Zenodo DOI, and optionally two or three
  notebooks reproducing Figs. 2, 5 and 8 with exercises.
- ORCID for the author: add `\orcid{...}` after the name in the `.tex` (needs the template's icon file) and give it again in the submission system.
- Cover letter: two paragraphs on the teaching gap the article fills and the reproducibility
  of every figure; mention that all numerical results were convergence-checked.
- Check the journal's current policy on anonymous review; if double-anonymous review is
  offered or required, remove the author line, addresses, acknowledgments and the repository
  URL from the reviewed copy.
- Article type: "Paper" (not "Letter").
