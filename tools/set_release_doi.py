#!/usr/bin/env python3
"""Write the version and Zenodo DOI of a release into every file that cites the code.

A release cannot contain its own version DOI: Zenodo mints it when the GitHub release is
created.  So a release is made in two steps:

  1. bump `version:` in CITATION.cff, commit, push, create the GitHub release vX.Y.Z;
     Zenodo archives it and shows the new version DOI;
  2. run this script with that DOI, rebuild the PDFs, commit and push:

         python3 tools/set_release_doi.py 1.2.0 10.5281/zenodo.NNNNNNNN

It rewrites
  CITATION.cff                          version, date-released, the list of release DOIs
  README.md                             the BibTeX entry and the sentence under it
  paper/qsol_cqed.tex                   reference `nguyen2026code`
  paper/ejp/jc_beyond_linear.tex        reference `nguyen2026code`
  paper/ejp/tools/make_aps_revtex.py    the same reference in Physical Review format
  paper/ejp/tools/make_ejp_iopjournal.py
and prints the commands that regenerate and rebuild the manuscripts.  The top-level `doi:`
of CITATION.cff stays the concept DOI, which always resolves to the latest release.
"""
import datetime
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONCEPT = '10.5281/zenodo.22806075'


def rewrite(path, rules):
    p = ROOT / path
    text = p.read_text()
    for pattern, replacement, expected in rules:
        text, n = re.subn(pattern, replacement, text)
        assert n == expected, f"{path}: {pattern!r} matched {n} times, expected {expected}"
    p.write_text(text)
    print(f"  updated {path}")


def main(version, doi, date):
    assert re.fullmatch(r"\d+\.\d+\.\d+", version), "version must look like 1.2.0"
    assert re.fullmatch(r"10\.5281/zenodo\.\d+", doi) and doi != CONCEPT, "give the VERSION DOI of the release"
    url = f"https://doi.org/{doi}"

    cff = (ROOT / 'CITATION.cff').read_text()
    extra = [] if doi in cff else [(r'(identifiers:\n(?:  - type: doi\n    value: \S+\n    description: "[^"]*"\n)+)',
                                    rf'\g<1>  - type: doi\n    value: {doi}\n    description: "Release v{version}"\n', 1)]
    rewrite('CITATION.cff', [(r"(?m)^version: \S+$", f"version: {version}", 1),
                             (r"(?m)^date-released: \S+$", f"date-released: {date}", 1)] + extra)
    rewrite('README.md', [
        (r"(  version = \{)[^}]*(\},)", rf"\g<1>{version}\g<2>", 1),
        (r"(  doi     = \{)[^}]*(\},)", rf"\g<1>{doi}\g<2>", 1),
        (r"(  url     = \{)https://doi\.org/10\.5281/zenodo\.\d+(\})", rf"\g<1>{url}\g<2>", 1),
        (r"The DOI above identifies release v[\d.]+, the version behind the manuscripts?'s? figures;",
         f"The DOI above identifies release v{version}, the version behind the manuscripts' figures;", 1),
    ])
    reference = [(r"(QSOL\\\\?_CQED\}?,? version )[\d.]+(,? (?:\(Zenodo\)|Zenodo \(\d{4}\),) \\\\?url\{)https://doi\.org/10\.5281/zenodo\.\d+(\})",
                  rf"\g<1>{version}\g<2>{url}\g<3>", 1)]
    for path in ('paper/qsol_cqed.tex', 'paper/ejp/jc_beyond_linear.tex',
                 'paper/ejp/tools/make_aps_revtex.py', 'paper/ejp/tools/make_ejp_iopjournal.py'):
        rewrite(path, reference)
    print(f"""
Release v{version}, {doi}.  Now regenerate and rebuild, then commit:

    python3 paper/ejp/tools/sync_supp_refs.py
    python3 paper/ejp/tools/make_aps_revtex.py
    (cd paper && tectonic qsol_cqed.tex)
    (cd paper/ejp && tectonic jc_beyond_linear.tex && tectonic jc_beyond_linear_supp.tex)
    (cd paper/aps && for f in jc_beyond_linear_aps jc_beyond_linear_aps_supp jc_beyond_linear_aps_research jc_beyond_linear_aps_research_supp; do tectonic $f.tex; done)
""")


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4) or '--help' in sys.argv:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else datetime.date.today().isoformat())
