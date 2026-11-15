"""Generate the APS (REVTeX 4.2) version of "Rotation and shear in the Jaynes-Cummings model".

Source (edit these, never the output):
    paper/ejp/jc_beyond_linear.tex        main text, IOP iopjournal class
    paper/ejp/jc_beyond_linear_supp.tex   supplementary material
Output (overwritten on every run):
    paper/aps/jc_beyond_linear_aps.tex        two-column `reprint' layout, as paper/qsol_cqed.tex
    paper/aps/jc_beyond_linear_aps_supp.tex   Supplemental Material, one column
    paper/aps/jc_main_refs_aps.tex            main-text numbers cited by the supplement
    paper/aps/jc_beyond_linear_aps_research*.tex, jc_main_refs_aps_research.tex
                                              the same three files as a research-journal cut: no teaching
                                              section, no exercises, no sentences addressed to students

    python3 paper/ejp/tools/make_aps_revtex.py
    cd paper/aps && tectonic jc_beyond_linear_aps.tex && tectonic jc_beyond_linear_aps_supp.tex

The physics, the numbers and the wording are those of the source.  What changes is the form:
  * class, title block, floats (figure*, ruledtabular), acknowledgments followed by a
    data-availability section (the order of the Physical Review journals);
  * cross-references in Physical Review style (Eq., Fig., Sec., Ref., Table with Roman numerals),
    the supplementary material cited as Supplemental Material through the reference list;
  * American spelling, serial commas in lists, double quotation marks, "et al." and "i.e.,";
  * displayed equations too wide for a 3.4 in column are set on two lines (EQUATIONS below);
    each keeps its single number, so the numbering is the same as in the source;
  * the reference list in Physical Review format (APS_REFS below), ordered by first citation.
Every substitution is checked: if the source changes in a place this script rewrites, the script
stops and says where, instead of producing a silently stale file.
"""
import hashlib
import pathlib
import re

EJP = pathlib.Path(__file__).resolve().parents[1]
OUT = EJP.parent / 'aps'
MAIN, SUPP = 'jc_beyond_linear', 'jc_beyond_linear_supp'

AUTHOR = r"""\author{Nguyen Khoi Nguyen}
\email{nknguyen@bu.edu}  % confirm the address to be published
\affiliation{Department of Electrical and Computer Engineering, Boston University, Boston, Massachusetts 02215, USA}
\affiliation{Department of Physics, Boston University, Boston, Massachusetts 02215, USA}
"""

MACROS = r"""\newcommand{\nbar}{\bar{n}}
\newcommand{\ket}[1]{|#1\rangle}
\newcommand{\bra}[1]{\langle #1|}
\newcommand{\braket}[2]{\langle #1|#2\rangle}
"""

# --------------------------------------------------------------------------------------------
# Reference list in Physical Review format.  Keys are those of the main text; the supplement
# uses the same entries under its own keys (SUPP_KEYS).
# --------------------------------------------------------------------------------------------
APS_REFS = {
    'jaynes1963': r"E.~T.~Jaynes and F.~W.~Cummings, Proc.\ IEEE \textbf{51}, 89 (1963).",
    'eberly1980': r"J.~H.~Eberly, N.~B.~Narozhny, and J.~J.~Sanchez-Mondragon, Phys.\ Rev.\ Lett.\ \textbf{44}, 1323 (1980).",
    'gerry2004': r"C.~C.~Gerry and P.~L.~Knight, \textit{Introductory Quantum Optics} (Cambridge University Press, Cambridge, 2004).",
    'haroche2006': r"S.~Haroche and J.-M.~Raimond, \textit{Exploring the Quantum: Atoms, Cavities, and Photons} (Oxford University Press, Oxford, 2006).",
    'gea-banacloche1990': r"J.~Gea-Banacloche, Phys.\ Rev.\ Lett.\ \textbf{65}, 3385 (1990).",
    'gea-banacloche1991': r"J.~Gea-Banacloche, Phys.\ Rev.\ A \textbf{44}, 5913 (1991).",
    'buzek1992': r"V.~Bu\v{z}ek, H.~Moya-Cessa, P.~L.~Knight, and S.~J.~D.~Phoenix, Phys.\ Rev.\ A \textbf{45}, 8190 (1992).",
    'phoenix1991': r"S.~J.~D.~Phoenix and P.~L.~Knight, Phys.\ Rev.\ A \textbf{44}, 6023 (1991).",
    'narozhny1981': r"N.~B.~Narozhny, J.~J.~Sanchez-Mondragon, and J.~H.~Eberly, Phys.\ Rev.\ A \textbf{23}, 236 (1981).",
    'fleischhauer1993': r"M.~Fleischhauer and W.~P.~Schleich, Phys.\ Rev.\ A \textbf{47}, 4258 (1993).",
    'feranchuk2009': r"I.~D.~Feranchuk and A.~V.~Leonov, Phys.\ Lett.\ A \textbf{373}, 517 (2009).",
    'pavlik2024': r"S.~I.~Pavlik, Phys.\ Scr.\ \textbf{99}, 015113 (2024).",
    'phoenix1988': r"S.~J.~D.~Phoenix and P.~L.~Knight, Ann.\ Phys.\ (N.Y.) \textbf{186}, 381 (1988).",
    'eiselt1991': r"J.~Eiselt and H.~Risken, Phys.\ Rev.\ A \textbf{43}, 346 (1991).",
    'larson2021': r"J.~Larson and T.~Mavrogordatos, \textit{The Jaynes--Cummings Model and Its Descendants} (IOP Publishing, Bristol, 2021).",
    'knight2024': r"P.~L.~Knight, C.~C.~Gerry, R.~J.~Birrittella, and P.~M.~Alsing, J.~Opt.\ Soc.\ Am.\ B \textbf{41}, C91 (2024).",
    'auffeves2003': r"A.~Auff\`{e}ves, P.~Maioli, T.~Meunier, S.~Gleyzes, G.~Nogues, M.~Brune, J.~M.~Raimond, and S.~Haroche, Phys.\ Rev.\ Lett.\ \textbf{91}, 230405 (2003).",
    'bild2023': r"M.~Bild \textit{et al.}, Science \textbf{380}, 274 (2023).",
    'meystre1982': r"P.~Meystre and M.~S.~Zubairy, Phys.\ Lett.\ A \textbf{89}, 390 (1982).",
    'kuklinski1988': r"J.~R.~Kukli\'{n}ski and J.~L.~Madajczyk, Phys.\ Rev.\ A \textbf{37}, 3175 (1988).",
    'woods1993': r"C.~W.~Woods and J.~Gea-Banacloche, J.~Mod.\ Opt.\ \textbf{40}, 2361 (1993).",
    'gea-banacloche1992': r"J.~Gea-Banacloche, Phys.\ Rev.\ A \textbf{46}, 7307 (1992).",
    'kennedy2010': r"R.~Kennedy, L.~Horstmeyer, A.~Dragan, and T.~Rudolph, Phys.\ Rev.\ A \textbf{82}, 054302 (2010).",
    'birrittella2015': r"R.~Birrittella, K.~Cheng, and C.~C.~Gerry, Opt.\ Commun.\ \textbf{354}, 286 (2015).",
    'assemat2019': r"F.~Assemat, D.~Grosso, A.~Signoles, A.~Facon, I.~Dotsenko, S.~Haroche, J.~M.~Raimond, M.~Brune, and S.~Gleyzes, Phys.\ Rev.\ Lett.\ \textbf{123}, 143605 (2019).",
    'gea-banacloche1993': r"J.~Gea-Banacloche, Phys.\ Rev.\ A \textbf{47}, 2221 (1993).",
    'meunier2006': r"T.~Meunier, A.~Le~Diffon, C.~Ruef, P.~Degiovanni, and J.-M.~Raimond, Phys.\ Rev.\ A \textbf{74}, 033802 (2006).",
    'munhoz2005': r"P.~P.~Munhoz and A.~Vidiella-Barranco, J.~Mod.\ Opt.\ \textbf{52}, 1557 (2005).",
    'kenfack2004': r"A.~Kenfack and K.~\.{Z}yczkowski, J.~Opt.~B: Quantum Semiclass.\ Opt.\ \textbf{6}, 396 (2004).",
    'johansson2012': r"J.~R.~Johansson, P.~D.~Nation, and F.~Nori, Comput.\ Phys.\ Commun.\ \textbf{183}, 1760 (2012).",
    'johansson2013': r"J.~R.~Johansson, P.~D.~Nation, and F.~Nori, Comput.\ Phys.\ Commun.\ \textbf{184}, 1234 (2013).",
    'lambert2026': r"N.~Lambert \textit{et al.}, Phys.\ Rep.\ \textbf{1153}, 1 (2026).",
    'walls1985': r"D.~F.~Walls and G.~J.~Milburn, Phys.\ Rev.\ A \textbf{31}, 2403 (1985).",
    'nguyen2026code': r"N.~K.~Nguyen, \textit{QSOL\_CQED}, version 1.2.0, Zenodo (2026), \url{https://doi.org/10.5281/zenodo.24262358}.",
}
SUPPLEMENT_REF = (r"See Supplemental Material at [URL will be inserted by publisher] for the derivations, the numerical "
                  r"methods and their measured errors, the limits of the asymptotic picture, other initial field states, "
                  r"maps of the cat state under cavity loss, and notes on the exercises.")
SUPP_KEYS = {'S-gea1991': 'gea-banacloche1991', 'S-gea1992': 'gea-banacloche1992', 'S-gea1993': 'gea-banacloche1993'}

# --------------------------------------------------------------------------------------------
# Displayed equations of the main text that do not fit a 3.4 in column: label -> (first ten
# hex digits of the SHA-1 of the source body, body for the two-column layout).
# --------------------------------------------------------------------------------------------
EQUATIONS = {
    'eq:exact': ('e5cc2b83a9', r"""\begin{aligned}
\ket{\psi(t)} &= \sum_{n=0}^{\infty} c_n\big[\cos\theta_n\,\ket{e,n} - i\sin\theta_n\,\ket{g,n+1}\big],\\
\theta_n &= gt\sqrt{n+1},\qquad c_n = e^{-\nbar/2}\frac{\alpha^n}{\sqrt{n!}}.
\end{aligned}"""),
    'eq:overlaps': ('0292146eb5', r"""\begin{aligned}
\langle\sigma_z\rangle &= \mathrm{Re}\,\braket{\Phi_-}{\Phi_+},\\
\langle(-1)^{\hat{n}}\rangle &= \mathrm{Re}\,\bra{\Phi_-}(-1)^{\hat{n}}\ket{\Phi_+}.
\end{aligned}"""),
    'eq:expansion': ('b207e89e69', r"""\begin{aligned}
\theta_n &\simeq gt\sqrt{\nbar}\,\Big[1 + \frac{\delta n}{2\nbar} - \frac{\delta n^2}{8\nbar^2} + \frac{\delta n^3}{16\nbar^3} - \dots\Big]\\
&= gt\sqrt{\nbar} + \varphi\,\delta n - \frac{\varphi\,\delta n^2}{4\nbar} + \frac{\varphi\,\delta n^3}{8\nbar^2} - \dots,
\end{aligned}"""),
    'eq:twobranch': ('860a5347bd', r"""\begin{aligned}
\ket{\psi(t)} &\simeq \frac{1}{\sqrt{2}}\Big[\ket{\Phi_-}\ket{\chi_-} + \ket{\Phi_+}\ket{\chi_+}\Big],\\
\ket{\chi_\mp} &= \frac{1}{\sqrt{2}}\left(\ket{e} \pm e^{\pm i\varphi}\ket{g}\right).
\end{aligned}"""),
    'eq:lemma': ('9386a7dcdf', r"""\begin{gathered}
\big\langle e^{-i\mu\,\delta n^2}\big\rangle = (1+2i\mu\nbar)^{-1/2}\\
\text{for $\delta n$ Gaussian with variance $\nbar$.}
\end{gathered}"""),
    'eq:revival': ('43f1b4dd96', r"""\begin{aligned}
&\big|\braket{\Phi_-}{\Phi_+}\big|_{t = kt_r} \simeq (1+\pi^2k^2)^{-1/4}\\
&\qquad = 0.551,\ 0.396,\ \dots\quad\text{for } k = 1,\ 2,\ \dots.
\end{aligned}"""),
    'eq:parity': ('2d6ddbd539', r"""\begin{aligned}
\langle(-1)^{\hat{n}}\rangle_{t_r/2} &\simeq -\,\frac{\cos(\pi\nbar-\eta)}{(1+\pi^2/4)^{1/4}},\\
\eta &= \tfrac{1}{2}\arctan\frac{\pi}{2} = 0.502,
\end{aligned}"""),
    'eq:nextorder': ('70e814f967', r"""\begin{aligned}
Z_k &= e^{2\pi ik(\nbar+1)}\,G(k)\left[1 + \frac{j(k)}{\nbar} + O(\nbar^{-2})\right],\\
j(k) &= -\frac{K(8K^2+3)}{96\,(K-i)^3},\qquad K = \pi k.
\end{aligned}"""),
    'eq:residual': ('0895cc9f64', r"""\begin{aligned}
\mathcal{P}(t_r/2) &= 1 - \frac{4+\pi^2}{32\,\nbar},\\
S(t_r/2) &\simeq h\!\left(\frac{4+\pi^2}{64\,\nbar}\right).
\end{aligned}"""),
}

# --------------------------------------------------------------------------------------------
# Wording.  (pattern, replacement, where) with where in {'main', 'supp', 'both'}; every regular
# expression must match at least once in the files it applies to.
# --------------------------------------------------------------------------------------------
WORDING = [
    # the supplementary material, cited through the reference list
    (r"\(supplementary section~(S\d)\)", r"(Sec.~\1 of the Supplemental Material)", 'main'),
    (r"\(supplementary sections~(S\d) and (S\d)\)", r"(Secs.~\1 and \2 of the Supplemental Material)", 'main'),
    (r"Supplementary section~(S\d)", r"Section~\1 of the Supplemental Material", 'main'),
    (r"supplementary section~(S\d)", r"Sec.~\1 of the Supplemental Material", 'main'),
    (r"the supplementary material", r"the Supplemental Material", 'main'),
    # cross-references (lower case only: a capital marks the start of a sentence and stays spelled out)
    (r"\bequation~\(", r"Eq.~(", 'both'),
    (r"\bequations \(B", r"Eqs.~(B", 'both'),
    (r"\bfigures~\\", r"Figs.~\\", 'supp'),
    (r"\bfigure~\\", r"Fig.~\\", 'both'),
    (r"\btable~\\", r"Table~\\", 'both'),
    (r"\bsections~(\\ref|S\d)", r"Secs.~\1", 'both'),
    (r"\bsection~(\\ref|S\d)", r"Sec.~\1", 'both'),
    (r"\breferences~\\cite", r"Refs.~\\cite", 'main'),
    (r"\breference~\\cite", r"Ref.~\\cite", 'both'),
    (r"\bappendix~B\b", r"Appendix~B", 'main'),
    # spelling and punctuation
    (r"\bcentrepiece\b", "centerpiece", 'main'),
    (r"\b([Cc])entred\b", r"\1entered", 'both'),
    (r"\b([Cc])entre(s?)\b", r"\1enter\2", 'both'),
    (r"\bgrey\b", "gray", 'both'),
    (r"\banalysed\b", "analyzed", 'both'),
    (r"\bneighbours\b", "neighbors", 'supp'),
    (r"\bper cent\b", "percent", 'main'),
    (r"\\textit\{et al\}", r"\\textit{et al.}", 'both'),
    (r"`phase spreading'", "``phase spreading''", 'main'),
    (r"i\.e\.\\ ", "i.e., ", 'supp'),
    # serial comma in lists of numbers: ", 25 and 100" -> ", 25, and 100"
    (r"(, \$?[+-]?\d[\d.]*(?:\\%)?\$?) and (\$?[+-]?\d)", r"\1, and \2", 'both'),
]

# Serial commas (Physical Review style): exact phrases, each of which must occur the stated
# number of times in the text it is applied to.  (phrase, replacement, where, count)
SERIAL = [
    (r"the disentanglement at $t_r/2$ and the cat state are due", r"the disentanglement at $t_r/2$, and the cat state are due", 'main', 1),
    (r"coherent states, density matrices and the partial trace", r"coherent states, density matrices, and the partial trace", 'main', 1),
    (r"in the Supplemental Material and in the figures", r"in the Supplemental Material, and in the figures", 'main', 1),
    (r"orthogonal, $\mathcal{P} = 1/2$ and $S = 1$~bit", r"orthogonal, $\mathcal{P} = 1/2$, and $S = 1$~bit", 'main', 1),
    (r"$t_r/4$, $t_r/2$ and $t_r$ (left to right)", r"$t_r/4$, $t_r/2$, and $t_r$ (left to right)", 'main', 1),
    (r"Cat fidelity, squeezing and negativity.", r"Cat fidelity, squeezing, and negativity.", 'main', 1),
    (r"the parity amplitude and the cat fidelity", r"the parity amplitude, and the cat fidelity", 'both', 1),
    (r"the skewness of the Poisson distribution and the difference", r"the skewness of the Poisson distribution, and the difference", 'main', 1),
    (r"$|Z_{1/2}|$, the parity and $F_{\rm cat}$", r"$|Z_{1/2}|$, the parity, and $F_{\rm cat}$", 'main', 1),
    (r"NumPy, SciPy, Matplotlib and QuTiP", r"NumPy, SciPy, Matplotlib, and QuTiP", 'main', 1),
    (r"$|G(\tfrac12)|$ and $|G(\tfrac14)|^2$", r"$|G(\tfrac12)|$, and $|G(\tfrac14)|^2$", 'main', 1),
    (r"5.9~dB at $\nbar = 100$ and 6.2~dB", r"5.9~dB at $\nbar = 100$, and 6.2~dB", 'main', 1),
    (r"Equation, figure and table numbers", r"Equation, figure, and table numbers", 'supp', 1),
    (r"$\xi = \delta n/\sqrt{\nbar}$ and $\epsilon", r"$\xi = \delta n/\sqrt{\nbar}$, and $\epsilon", 'supp', 1),
    (r"$4\times10^{-7}$ at 30 and $4\times10^{-8}$ at 34", r"$4\times10^{-7}$ at 30, and $4\times10^{-8}$ at 34", 'supp', 1),
    (r"a slow, smooth and entirely spurious", r"a slow, smooth, and entirely spurious", 'supp', 1),
    (r"80~s at 64 and 340~s at 100", r"80~s at 64, and 340~s at 100", 'supp', 1),
    (r"thermal, squeezed vacuum and Fock)", r"thermal, squeezed vacuum, and Fock)", 'supp', 1),
    (r"The rotation, the shear and the squeezing", r"The rotation, the shear, and the squeezing", 'supp', 1),
    (r"(e), (f) squeezed-vacuum and (g), (h) Fock", r"(e), (f) squeezed-vacuum, and (g), (h) Fock", 'supp', 1),
]


def sha(s):
    return hashlib.sha1(s.encode()).hexdigest()[:10]


def strip_comments(t):
    return re.sub(r"(?m)^[ \t]*%.*\n", "", re.sub(r"(?<!\\)[ \t]+%.*", "", t))


def take(pattern, text, what):
    m = re.search(pattern, text, re.S)
    assert m, f"could not find {what} in the source"
    return m


def braced(text, macro):
    """Argument of \\macro{...} (balanced braces)."""
    i = text.index('\\' + macro + '{') + len(macro) + 2
    depth, j = 1, i
    while depth:
        depth += {'{': 1, '}': -1}.get(text[j], 0)
        j += 1
    return text[i:j - 1]


def wording(t, where):
    for pat, rep, scope in WORDING:
        if scope not in (where, 'both'):
            continue
        t, n = re.subn(pat, rep, t)
        assert n > 0, f"wording rule {pat!r} no longer matches the {where} text: update WORDING"
    for phrase, rep, scope, count in SERIAL:
        if scope not in (where, 'both'):
            continue
        assert t.count(phrase) == count, f"serial-comma phrase {phrase!r}: {t.count(phrase)} occurrences in {where}, expected {count}"
        t = t.replace(phrase, rep)
    for stale in ('upplementary', 'equation~(', 'figure~\\', 'reference~\\cite', 'et al}', 'centre', 'grey'):
        hits = [m.start() for m in re.finditer(re.escape(stale), t) if not (stale == 'equation~(' and t[m.start() - 1] == 'E')
                and not (stale == 'figure~\\' and t[m.start() - 1] == 'F')]
        assert not hits, f"{stale!r} survived in the {where} text near: {t[max(0, hits[0] - 60):hits[0] + 40]!r}"
    return t


def floats(t, two_column):
    star = '*' if two_column else ''

    def table(m):
        cap, spec, rows, label = m.group(1), m.group(2), m.group(3).strip(), m.group(4)
        assert rows.startswith('\\hline') and rows.endswith('\\hline'), label
        rows = rows[len('\\hline'):-len('\\hline')].strip()
        return (f"\\begin{{table{star}}}[t]\n\\caption{{{cap}}}\n\\label{{{label}}}\n\\begin{{ruledtabular}}\n"
                f"\\begin{{tabular}}{{{spec}}}\n{rows}\n\\end{{tabular}}\n\\end{{ruledtabular}}\n\\end{{table{star}}}")
    t, n_tab = re.subn(r"\\begin\{table\}\n\\centering\n\\caption\{(.*?)\}\n\\begin\{tabular\}\{([^}]*)\}\n(.*?)\n\\end\{tabular\}\n"
                       r"\\label\{([^}]*)\}\n\\end\{table\}", table, t, flags=re.S)

    def figure(m):
        width = {'\\textwidth': '6in', '0.92\\textwidth': '5.5in'}[m.group(1)]
        return (f"\\begin{{figure{star}}}[t]\n\\centering\n\\includegraphics[width={width}]{{{m.group(2)}}}\n"
                f"\\caption{{{m.group(3)}}}\n\\label{{{m.group(4)}}}\n\\end{{figure{star}}}")
    t, n_fig = re.subn(r"\\begin\{figure\}\n\\centering\n\\includegraphics\[width=([^\]]*)\]\{([^}]*)\}\n\\caption\{(.*?)\}\n"
                       r"\\label\{([^}]*)\}\n\\end\{figure\}", figure, t, flags=re.S)
    assert '\\begin{table}\n\\centering' not in t and '\\begin{figure}\n\\centering' not in t and '\\hline\n\\end{tabular}' not in t
    return t, n_tab, n_fig


def equations(t):
    seen = set()

    def one(m):
        body, label = m.group(1), m.group(2)
        if label not in EQUATIONS:
            return m.group(0)
        want, new = EQUATIONS[label]
        assert sha(body) == want, (f"equation {label} changed in the source (hash {sha(body)}, expected {want!r}): "
                                   "bring its two-line form in EQUATIONS up to date, then update the hash")
        seen.add(label)
        return f"\\begin{{equation}}\n{new}\n\\label{{{label}}}\n\\end{{equation}}"
    t = re.sub(r"\\begin\{equation\}\n(.*?)\n\\label\{([^}]*)\}\n\\end\{equation\}", one, t, flags=re.S)
    assert seen == set(EQUATIONS), set(EQUATIONS) - seen
    return t


def check_reference(key, iop, aps):
    """The Physical Review entry must carry the year, volume, page and surnames of the source entry."""
    iop = iop.replace('{\\it et al}', '')
    year = take(r"\b((?:19|20)\d\d)\b", iop, f"year of {key}").group(1)
    assert f"{year})" in aps or f"({year})" in aps or f", {year})" in aps, (key, year, aps)
    vol = re.search(r"\{\\bf (\w+)\}\s+(\w+)", iop)
    if vol:
        assert f"\\textbf{{{vol.group(1)}}}, {vol.group(2)} " in aps, (key, vol.groups(), aps)
    authors = iop[:iop.index(year)]
    for person in re.split(r", | and ", authors.strip()):
        m = re.match(r"^(.*?)\s+(?:[A-Z](?:-[A-Z])?)(?: [A-Z])*$", person.strip())
        if m:
            assert m.group(1) in aps.replace('~', ' '), (key, m.group(1), aps)


def bibliography(body, source_bib, extra=None, keymap=None):
    """thebibliography in order of first citation; the cited keys must be exactly the available ones."""
    source = dict(re.findall(r"\\bibitem\{([^}]*)\}\s*(.*)", source_bib))
    entries = {}
    for key, iop in source.items():
        aps = APS_REFS[(keymap or {}).get(key, key.removeprefix('S-'))]
        check_reference(key, iop, aps)
        entries[key] = aps
    entries.update(extra or {})
    cited = []
    for m in re.finditer(r"\\cite\{([^}]*)\}", body):
        for key in (k.strip() for k in m.group(1).split(',')):
            if key not in cited:
                cited.append(key)
    assert set(cited) == set(entries), set(cited) ^ set(entries)
    return ("\\begin{thebibliography}{99}\n\n" + "\n\n".join(f"\\bibitem{{{k}}}\n{entries[k]}" for k in cited)
            + "\n\n\\end{thebibliography}\n"), len(cited)


def roman(n):
    out = ''
    for value, sym in ((10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')):
        while n >= value:
            out, n = out + sym, n - value
    return out


def main_refs(body):
    """Numbers of the main text's equations, figures, tables and sections, for \\mainref in the supplement.
    Tables and sections carry Roman numerals in the APS layout."""
    refs = {}
    for env in ('equation', 'figure', 'table'):
        blocks = re.findall(r"\\begin\{%s\*?\}(.*?)\\end\{%s\*?\}" % (env, env), body, re.S)
        for number, blk in enumerate(blocks, 1):
            labels = re.findall(r"\\label\{([^}]*)\}", blk)
            assert len(labels) == 1, (env, number, labels)
            refs[labels[0]] = roman(number) if env == 'table' else number
    for number, m in enumerate(re.finditer(r"\\section\{.*?\}\s*\\label\{([^}]*)\}", body, re.S), 1):
        refs[m.group(1)] = roman(number)
    lines = ["% Generated by tools/make_aps_revtex.py from jc_beyond_linear.tex. Do not edit.", "\\makeatletter",
             "\\newcommand{\\mainref}[1]{\\@ifundefined{mainref@#1}{\\textbf{??}}{\\@nameuse{mainref@#1}}}"]
    lines += [f"\\@namedef{{mainref@{k}}}{{{v}}}" for k, v in refs.items()] + ["\\makeatother", ""]
    return '\n'.join(lines), refs


HEADER = """%% ==========================================================================
%% GENERATED by paper/ejp/tools/make_aps_revtex.py from paper/ejp/{src}.tex.
%% Do not edit this file: edit the source and run the script again.{variant}
%% Build: tectonic {out}.tex
%% ==========================================================================
"""


RESEARCH_NOTE = """
%% RESEARCH CUT: without the teaching section, the exercises and the sentences addressed to
%% students (CUT_* in the script)."""


# --------------------------------------------------------------------------------------------
# Research cut: the same article without the teaching section, the exercises and the sentences
# addressed to students.  Applied to the APS text; (phrase, replacement), each phrase must occur
# exactly once.  The section "Using this material in teaching" and the supplement's "Notes on
# the exercises" are removed whole.
# --------------------------------------------------------------------------------------------
CUT_KEYWORDS = [(r", computational physics education", "")]
CUT_ABSTRACT = [(r" The paper is aimed at advanced undergraduate and beginning graduate students; exercises and a "
                 r"Python script for every figure are provided.", "")]
CUT_MAIN = [
    (r"is where most students first meet the quantized field doing something a classical field cannot.",
     r"is the simplest setting in which the quantized field does something a classical field cannot."),
    (r"A student who then runs a simulation sees a good deal that this argument does not explain.",
     r"A simulation shows a good deal that this argument does not explain."),
    (r" A student can test each of these to three or four digits.",
     r" Each of these can be tested against the exact solution to three or four digits."),
    (r"We assume familiarity with second quantization of one mode, coherent states, density matrices, and the partial trace. ", ""),
    (r", Sec.~\ref{sec:loss} adds cavity loss, and Sec.~\ref{sec:teaching} suggests how to use the material in a course.",
     r", and Sec.~\ref{sec:loss} adds cavity loss."),
    (r"and a short calculation (exercise 2) gives", r"and a short calculation gives"),
    ("\n\nFor teaching, the value of the exercise is that it closes a loop that a simulation alone leaves open: the student "
     "predicts a number, computes it, and has to account for the difference.", ""),
]
CUT_SUPP = [
    (r" Section~S6 gives notes on the exercises.", ""),
    (r"~\cite{S-gea1991} (exercise 7).", r"~\cite{S-gea1991}."),
]
SUPPLEMENT_REF_RESEARCH = (r"See Supplemental Material at [URL will be inserted by publisher] for the derivations, the "
                           r"numerical methods and their measured errors, the limits of the asymptotic picture, other "
                           r"initial field states, and maps of the cat state under cavity loss.")
TEACHING_WORDS = r"(?i)student|teach|exercise|lecture|undergraduate|classroom|pedagog"


def cut(text, cuts, what):
    for phrase, rep in cuts:
        assert text.count(phrase) == 1, f"research cut: {phrase[:70]!r} occurs {text.count(phrase)} times in the {what}"
        text = text.replace(phrase, rep)
    return text


def build_main(research=False):
    tag = '_aps_research' if research else '_aps'
    src = strip_comments((EJP / f'{MAIN}.tex').read_text())
    doc = src.split('\\begin{document}', 1)[1]
    title = braced(doc, 'title')
    assert title.count(': what the') == 1 and title.endswith('cat states and entanglement'), title
    title = title.replace(': what the', ': What the').replace('cat states and entanglement', 'cat states, and entanglement')
    abstract = take(r"\\begin\{abstract\}\n(.*?)\n\\end\{abstract\}", doc, 'abstract').group(1)
    keywords = braced(doc, 'keywords')
    body = take(r"(\\section\{Introduction\}.*?)\n\\ack\{", doc, 'body').group(1)
    ack, data = braced(doc, 'ack'), braced(doc, 'data')
    source_bib = doc[doc.index('\\begin{thebibliography}'):]

    body, n_tab, n_fig = floats(equations(body), two_column=True)
    text = wording('\n'.join([abstract, '%%SPLIT%%', body, '%%SPLIT%%', data, '%%SPLIT%%', ack]), 'main')
    abstract, body, data, ack = (part.strip('\n') for part in text.split('%%SPLIT%%'))
    if research:
        keywords, abstract, body = cut(keywords, CUT_KEYWORDS, 'keywords'), cut(abstract, CUT_ABSTRACT, 'abstract'), cut(body, CUT_MAIN, 'main text')
        i, j = body.index('\\section{Using this material in teaching}'), body.index('\\section{Conclusions}')
        assert 0 < i < j and '\\begin{equation}' not in body[i:j] and '\\begin{figure' not in body[i:j] and '\\begin{table' not in body[i:j]
        body = body[:i] + body[j:]
        left = re.findall(TEACHING_WORDS, abstract + body + keywords)
        assert not left, f"teaching vocabulary left in the research cut: {left}"
    assert 'Supplemental Material' not in abstract
    assert body.count('the Supplemental Material treat') == 1
    body = body.replace('the Supplemental Material treat', 'the Supplemental Material~\\cite{supplement} treat')
    bib, n_ref = bibliography(body + data, source_bib, extra={'supplement': SUPPLEMENT_REF_RESEARCH if research else SUPPLEMENT_REF})

    out = HEADER.format(src=MAIN, out=MAIN + tag, variant=RESEARCH_NOTE if research else '') + r"""\documentclass[%
 reprint,
 amsmath,amssymb,
 aps,
 floatfix,
]{revtex4-2}

\usepackage{graphicx}
\graphicspath{{../figures/}{../../figures/}}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
""" + MACROS + r"""
\begin{document}

\title{""" + title + "}\n\n" + AUTHOR + r"""
\date{\today}

\begin{abstract}
""" + abstract + r"""
\end{abstract}

\keywords{""" + keywords + r"""}

\maketitle

""" + body + r"""

% Keep the two closing sections together: a short block left alone on the last page is set
% in the right-hand column only when REVTeX balances the columns above the reference list.
\vspace{0pt plus 2.6in}\pagebreak[3]\vspace{0pt plus -2.6in}
\begin{acknowledgments}
""" + ack + r"""
\end{acknowledgments}

\section*{Data availability}

""" + data + r"""

""" + bib + "\n\\end{document}\n"
    (OUT / f'{MAIN}{tag}.tex').write_text(out)
    refs_tex, refs = main_refs(body)
    (OUT / f'jc_main_refs{tag}.tex').write_text(refs_tex)
    print(f"{MAIN}{tag}.tex: {n_fig} figures, {n_tab} tables, {len(EQUATIONS)} equations reset on two lines, {n_ref} references")
    return title, refs


def build_supp(title, refs, research=False):
    tag = '_aps_research' if research else '_aps'
    src = strip_comments((EJP / f'{SUPP}.tex').read_text())
    doc = src.split('\\begin{document}', 1)[1]
    intro = take(r"\\begin\{abstract\}\n(.*?)\n\\end\{abstract\}", doc, 'abstract').group(1)
    body = take(r"(\\section\{Second order: derivations\}.*?)\n\\begin\{thebibliography\}", doc, 'body').group(1)
    source_bib = doc[doc.index('\\begin{thebibliography}'):]

    body, n_tab, n_fig = floats(body, two_column=False)
    text = wording(intro + '\n%%SPLIT%%\n' + body, 'supp')
    intro, body = (part.strip('\n') for part in text.split('%%SPLIT%%'))
    if research:
        intro, body = cut(intro, CUT_SUPP[:1], 'supplement'), cut(body, CUT_SUPP[1:], 'supplement')
        i = body.index('\\clearpage\n\\section{Notes on the exercises}')
        assert body.count('\\section{') == 6 and body[i:].count('\\section{') == 1, "the notes on the exercises must be the last section"
        body = body[:i].rstrip('\n')
        left = re.findall(TEACHING_WORDS, intro + body)
        assert not left, f"teaching vocabulary left in the research supplement: {left}"
    used = set(re.findall(r"\\mainref\{([^}]*)\}", intro + body))
    assert used <= set(refs), used - set(refs)
    bib, n_ref = bibliography(body, source_bib, keymap=SUPP_KEYS)

    out = HEADER.format(src=SUPP, out=MAIN + tag + '_supp', variant=RESEARCH_NOTE if research else '') + r"""\documentclass[%
 reprint,onecolumn,
 amsmath,amssymb,
 aps,
 floatfix,
]{revtex4-2}

\usepackage{graphicx}
\graphicspath{{../figures/}{../../figures/}}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
""" + MACROS + r"\input{jc_main_refs" + tag + "}\n" + r"""\setlength{\emergencystretch}{2em}
\raggedbottom
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thesubsection}{\thesection.\arabic{subsection}}
\renewcommand{\theequation}{S\arabic{equation}}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}

\begin{document}

\title{Supplemental Material for ``""" + title + "''}\n\n" + AUTHOR + r"""
\date{\today}

\maketitle

""" + intro + "\n\n" + body + "\n\n" + bib + "\n\\end{document}\n"
    (OUT / f'{MAIN}{tag}_supp.tex').write_text(out)
    print(f"{MAIN}{tag}_supp.tex: {n_fig} figures, {n_tab} tables, {n_ref} references, {len(used)} main-text labels cited")


if __name__ == '__main__':
    import sys
    if '--hashes' in sys.argv:                     # print the hashes of the source equations listed in EQUATIONS
        doc = strip_comments((EJP / f'{MAIN}.tex').read_text())
        for body, label in re.findall(r"\\begin\{equation\}\n(.*?)\n\\label\{([^}]*)\}\n\\end\{equation\}", doc, re.S):
            if label in EQUATIONS:
                print(f"    {label!r}: {sha(body)!r}")
        raise SystemExit
    OUT.mkdir(exist_ok=True)
    for research in (False, True):
        build_supp(*build_main(research), research)
