"""Convert paper/qsol_cqed.tex (REVTeX 4.2) into IOP Publishing's current
iopjournal.cls template (2024) for the European Journal of Physics.

Content is preserved verbatim; only markup changes:
  * iopjournal title block: \\articletype, \\title, \\author, \\affil, \\email,
    \\keywords, abstract (no \\maketitle, no short title, no \\submitto)
  * Figures 1-8 placed after the paragraph that first cites them; Figs. 9-12
    (uncited in the text) omitted by request
  * tables in the template's plain \\hline style
  * \\bm -> \\boldsymbol (amsmath); \\ack{...} and \\data{...} back matter
  * bibliography in IOP numeric style (thebibliography prints its own heading)
"""
import re, pathlib, sys
from PIL import Image

ROOT = pathlib.Path('/Users/alanknguyen/Documents/qsol_cqed')
SRC = (ROOT / 'paper' / 'qsol_cqed.tex').read_text()
OUTDIR = ROOT / 'paper' / 'ejp'; OUTDIR.mkdir(exist_ok=True)
OUT = OUTDIR / 'qsol_cqed_ejp.tex'
OMIT = {'fig:rabi_collapse_revival', 'fig:wigner_quantum_states', 'fig:wigner_2d_quantum_states', 'fig:mollow_triplet'}

# ---------------------------------------------------------------- pieces of the source
abstract = re.search(r'\\begin\{abstract\}\n(.*?)\n\\end\{abstract\}', SRC, re.S).group(1).strip()
keywords = re.search(r'\\keywords\{((?:[^{}]|\{[^{}]*\})*)\}', SRC).group(1)
body = SRC.split('\\maketitle', 1)[1].split('\\newpage', 1)[0]
figpart = SRC.split('\\newpage', 1)[1].split('\\begin{acknowledgments}', 1)[0]
ack = re.search(r'\\begin\{acknowledgments\}\n(.*?)\n\\end\{acknowledgments\}', SRC, re.S).group(1).strip()
bibkeys = re.findall(r'\\bibitem\{([^}]*)\}', SRC)

# ---------------------------------------------------------------- figures
figures = []
for blk in re.findall(r'\\begin\{figure\*?\}\[t\](.*?)\\end\{figure\*?\}', figpart, re.S):
    fname = re.search(r'\\includegraphics\[[^\]]*\]\{([^}]*)\}', blk).group(1)
    cap = re.search(r'\\caption\{(.*)\}\s*\n\s*\\label\{([^}]*)\}', blk, re.S)
    caption, label = cap.group(1).strip(), cap.group(2)
    w, h = Image.open(ROOT / 'paper' / 'figures' / fname).size
    width = r'0.62\textwidth' if h / w > 1.2 else (r'0.85\textwidth' if h / w > 0.75 else r'\textwidth')
    figures.append(dict(file=fname, caption=caption, label=label, width=width))
assert len(figures) == 12

def fig_env(f):
    return ('\\begin{figure}\n\\centering\n'
            f'\\includegraphics[width={f["width"]}]{{{f["file"]}}}\n'
            f'\\caption{{{f["caption"]}}}\n'
            f'\\label{{{f["label"]}}}\n'
            '\\end{figure}')

# ---------------------------------------------------------------- body
body = re.sub(r'%% =+\n%% [^\n]*\n%% =+\n', '', body)
body = body.replace('\\bm{', '\\boldsymbol{')

def convert_table(m):
    blk = m.group(0)
    blk = blk.replace('\\begin{table}[h]', '\\begin{table}')
    blk = blk.replace('\\hline\\hline', '\\hline')
    # template style: caption above, label after the tabular
    lab = re.search(r'\n\\label\{([^}]*)\}\n', blk).group(1)
    blk = blk.replace(f'\n\\label{{{lab}}}\n', '\n', 1)
    blk = blk.replace('\\end{tabular}\n', f'\\end{{tabular}}\n\\label{{{lab}}}\n', 1)
    return blk

body = re.sub(r'\\begin\{table\}\[h\].*?\\end\{table\}', convert_table, body, flags=re.S)

paras = body.split('\n\n'); placed = set()
for f in figures:
    if f['label'] in OMIT: continue
    for i, p in enumerate(paras):
        if f'\\ref{{{f["label"]}}}' in p and not p.lstrip().startswith('\\begin{figure}'):
            paras[i] = p + '\n\n' + fig_env(f); placed.add(f['label']); break
assert placed == {f['label'] for f in figures} - OMIT, placed
body = '\n\n'.join(paras).strip('\n')

# ---------------------------------------------------------------- bibliography (IOP numeric style)
IOP = {
 'rabi1937':   r"Rabi I I 1937 {\it Phys. Rev.} {\bf 51} 652",
 'scully1997': r"Scully M O and Zubairy M S 1997 {\it Quantum Optics} (Cambridge: Cambridge University Press)",
 'allen1975':  r"Allen L and Eberly J H 1975 {\it Optical Resonance and Two-Level Atoms} (New York: Dover)",
 'gerry2004':  r"Gerry C C and Knight P L 2004 {\it Introductory Quantum Optics} (Cambridge: Cambridge University Press)",
 'jaynes1963': r"Jaynes E T and Cummings F W 1963 {\it Proc. IEEE} {\bf 51} 89",
 'eberly1980': r"Eberly J H, Narozhny N B and Sanchez-Mondragon J J 1980 {\it Phys. Rev. Lett.} {\bf 44} 1323",
 'brune1996':  r"Brune M, Hagley E, Dreyer J, Ma\^{\i}tre X, Maali A, Wunderlich C, Raimond J M and Haroche S 1996 {\it Phys. Rev. Lett.} {\bf 77} 4887",
 'rempe1987':  r"Rempe G, Walther H and Klein N 1987 {\it Phys. Rev. Lett.} {\bf 58} 353",
 'haroche2006': r"Haroche S and Raimond J-M 2006 {\it Exploring the Quantum: Atoms, Cavities, and Photons} (Oxford: Oxford University Press)",
 'blais2021':  r"Blais A, Grimsmo A L, Girvin S M and Wallraff A 2021 {\it Rev. Mod. Phys.} {\bf 93} 025005",
 'glauber1963': r"Glauber R J 1963 {\it Phys. Rev.} {\bf 131} 2766",
 'sudarshan1963': r"Sudarshan E C G 1963 {\it Phys. Rev. Lett.} {\bf 10} 277",
 'aasi2013':   r"Aasi J {\it et al} (LIGO Scientific Collaboration) 2013 {\it Nat. Photon.} {\bf 7} 613",
 'tse2019':    r"Tse M {\it et al} (LIGO Scientific Collaboration) 2019 {\it Phys. Rev. Lett.} {\bf 123} 231107",
 'wigner1932': r"Wigner E 1932 {\it Phys. Rev.} {\bf 40} 749",
 'hudson1974': r"Hudson R L 1974 {\it Rep. Math. Phys.} {\bf 6} 249",
 'kenfack2004': r"Kenfack A and \.{Z}yczkowski K 2004 {\it J. Opt. B: Quantum Semiclass. Opt.} {\bf 6} 396",
 'mollow1969': r"Mollow B R 1969 {\it Phys. Rev.} {\bf 188} 1969",
 'thompson1992': r"Thompson R J, Rempe G and Kimble H J 1992 {\it Phys. Rev. Lett.} {\bf 68} 1132",
 'lindblad1976': r"Lindblad G 1976 {\it Commun. Math. Phys.} {\bf 48} 119",
 'wiseman2009': r"Wiseman H M and Milburn G J 2009 {\it Quantum Measurement and Control} (Cambridge: Cambridge University Press)",
 'zurek2003':  r"Zurek W H 2003 {\it Rev. Mod. Phys.} {\bf 75} 715",
 'nielsen2000': r"Nielsen M A and Chuang I L 2000 {\it Quantum Computation and Quantum Information} (Cambridge: Cambridge University Press)",
 'phoenix1988': r"Phoenix S J D and Knight P L 1988 {\it Ann. Phys. (N.Y.)} {\bf 186} 381",
 'johansson2012': r"Johansson J R, Nation P D and Nori F 2012 {\it Comput. Phys. Commun.} {\bf 183} 1760",
 'johansson2013': r"Johansson J R, Nation P D and Nori F 2013 {\it Comput. Phys. Commun.} {\bf 184} 1234",
 'gea-banacloche1990': r"Gea-Banacloche J 1990 {\it Phys. Rev. Lett.} {\bf 65} 3385",
 'gea-banacloche1991': r"Gea-Banacloche J 1991 {\it Phys. Rev.} A {\bf 44} 5913",
 'phoenix1991': r"Phoenix S J D and Knight P L 1991 {\it Phys. Rev.} A {\bf 44} 6023",
 'vidal2002':  r"Vidal G and Werner R F 2002 {\it Phys. Rev.} A {\bf 65} 032314",
 'plenio2005': r"Plenio M B 2005 {\it Phys. Rev. Lett.} {\bf 95} 090503",
 'kimble1977': r"Kimble H J, Dagenais M and Mandel L 1977 {\it Phys. Rev. Lett.} {\bf 39} 691",
 'ekert1991':  r"Ekert A K 1991 {\it Phys. Rev. Lett.} {\bf 67} 661",
 'hofheinz2009': r"Hofheinz M, Wang H, Ansmann M, Bialczak R C, Lucero E, Neeley M, O'Connell A D, Sank D, Wenner J, Martinis J M and Cleland A N 2009 {\it Nature} {\bf 459} 546",
 'ofek2016':   r"Ofek N {\it et al} 2016 {\it Nature} {\bf 536} 441",
 'meekhof1996': r"Meekhof D M, Monroe C, King B E, Itano W M and Wineland D J 1996 {\it Phys. Rev. Lett.} {\bf 76} 1796",
 'walls2008':  r"Walls D F and Milburn G J 2008 {\it Quantum Optics} 2nd edn (Berlin: Springer)",
}
assert set(bibkeys) == set(IOP)
bib = '\n'.join(f'\\bibitem{{{k}}} {IOP[k]}' for k in bibkeys)

title = ("Quantum states of light in cavity QED: a computational study of Wigner function "
         "dynamics and atom--field entanglement in the Jaynes--Cummings model")

doc = f"""% ==========================================================================
% European Journal of Physics (IOP Publishing) version of qsol_cqed.tex,
% using IOP's current iopjournal.cls template (2024).  On Overleaf: create the
% project from the "IOP Publishing" journal template (it ships iopjournal.cls and
% the orcid icon), replace main.tex with this file, and upload the figures into
% a folder named figures/.  Content is identical to the REVTeX version.
% Generated by scratchpad/make_ejp_iopjournal.py from paper/qsol_cqed.tex.
% Options: \\documentclass[anonymous]{{iopjournal}} strips author details,
% acknowledgments and the data statement for double-anonymous review.
% ==========================================================================
\\documentclass{{iopjournal}}
% graphicx, fancyhdr, xcolor and hyperref are loaded by the class.
\\usepackage{{amsmath}}
\\usepackage{{amssymb}}
\\graphicspath{{{{figures/}}{{../figures/}}}}
% running head placeholders (volume/year/article number are set by the journal)
\\fancyhead[L]{{{{\\small \\sf IOP Publishing}}\\hspace{{5mm}} {{\\it Eur. J. Phys.}} {{\\bf vv}} (yyyy) aaaaaa}}
\\fancyhead[R]{{N K Nguyen}}

\\begin{{document}}

\\articletype{{Paper}}

\\title{{{title}}}

% To add an ORCID iD use \\orcid{{0000-0000-0000-0000}} after the name; it needs the
% "orcid" icon file that ships with the Overleaf IOP template.
\\author{{Nguyen Khoi Nguyen$^{{1,2,*}}$}}

\\affil{{$^1$Department of Electrical and Computer Engineering, Boston University, Boston, MA 02215, United States of America}}

\\affil{{$^2$Physics Department, Boston University, Boston, MA 02215, United States of America}}

\\affil{{$^*$Author to whom any correspondence should be addressed.}}

\\email{{nknguyen@bu.edu}}  % confirm the address to be published

\\keywords{{{keywords}}}

\\begin{{abstract}}
{abstract}
\\end{{abstract}}

{body}

\\ack{{{ack}}}

% \\funding{{This work received no external funding.}}
% \\roles{{N K N carried out all the work reported.}}

\\data{{The simulation code that generates every figure in this article, together with the raw sweep data, is openly available at \\url{{https://github.com/alanknguyen/QSOL_CQED}}. All results were obtained with QuTiP 5.2.3.}}

% \\suppdata{{...}}

% thebibliography prints the "References" heading itself in this (article-based) class
\\begin{{thebibliography}}{{99}}
{bib}
\\end{{thebibliography}}

\\end{{document}}
"""
OUT.write_text(doc)

# ---------------------------------------------------------------- static checks
envs = re.findall(r'\\begin\{([a-zA-Z*]+)\}', doc); ends = re.findall(r'\\end\{([a-zA-Z*]+)\}', doc)
assert all(envs.count(k) == ends.count(k) for k in set(envs) | set(ends))
labels = set(re.findall(r'\\label\{([^}]*)\}', doc)); refs = set(re.findall(r'\\ref\{([^}]*)\}', doc))
assert refs <= labels, refs - labels
cites = set(k.strip() for m in re.findall(r'\\cite\{([^}]*)\}', doc) for k in m.split(','))
assert cites <= set(bibkeys)
t = re.sub(r'(?<!\\)%.*', '', doc).replace('\\{', '').replace('\\}', '')
assert t.count('{') == t.count('}'), (t.count('{'), t.count('}'))
abs_words = len(re.sub(r'\\[a-zA-Z]+|[{}$^_]', ' ', abstract).split())
print(f"wrote {OUT} ({len(doc)} chars); {len(placed)} figures placed, {len(OMIT)} omitted; "
      f"{len(refs)} refs and {len(cites)} cite keys resolved; abstract ~{abs_words} words")
