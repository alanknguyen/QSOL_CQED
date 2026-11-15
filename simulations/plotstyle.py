# plotstyle.py
# Shared Matplotlib style and small plotting helpers for every figure and
# animation in the repository.
#
# Figures are drawn at their final printed size (IOP single-column text width
# is 153 mm = 6.0 in), so the font sizes below are the sizes on the page:
# 8-10 pt text, panels labelled (a), (b), ... and no figure-level titles (the
# caption carries that information).  Line art is saved as vector PDF with
# embedded TrueType fonts; density maps are rasterised inside the PDF so the
# files stay small.
#
# Nguyen Khoi Nguyen (Alan), Boston University

import io
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PAPER_FIG = ROOT / 'paper' / 'figures'
FIG_DIR = ROOT / 'figures'
ANIM_DIR = ROOT / 'animations'
SUPP_DIR = ROOT / 'paper' / 'supplementary'
for _d in (PAPER_FIG, FIG_DIR, ANIM_DIR):
    _d.mkdir(parents=True, exist_ok=True)

FULL_W = 6.0      # inches: full text width of the IOP layout (153 mm)
HALF_W = 3.4      # inches: one column of a two-column layout
W_MAX = 1.0 / np.pi   # largest |W| any state can reach in the a = (x + ip)/sqrt(2) convention

# Colours chosen to stay distinguishable in greyscale and for colour-blind
# readers; every multi-line plot also varies the line style.
COLORS = {
    'blue':   '#2166ac',
    'red':    '#b2182b',
    'green':  '#1b7837',
    'orange': '#e08214',
    'purple': '#762a83',
    'grey':   '#6e6e6e',
}
GUIDE_HALF = dict(color='#b2182b', lw=0.7, ls=':')    # t_r / 2
GUIDE_FULL = dict(color='#4d4d4d', lw=0.7, ls='--')   # t_r


def use_style():
    """Apply the repository-wide rcParams.  Call once at the top of a script."""
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans'],
        'mathtext.fontset': 'dejavusans',
        'font.size': 8.5,
        'axes.labelsize': 9.5,
        'axes.titlesize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'legend.frameon': False,
        'legend.handlelength': 2.2,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'lines.linewidth': 1.0,
        'lines.markersize': 4,
        'figure.dpi': 110,
        'savefig.dpi': 300,
        'pdf.fonttype': 42,       # embed TrueType so text stays text
        'ps.fonttype': 42,
        'figure.constrained_layout.use': True,
    })


def sequential_colors(n, cmap='viridis', lo=0.0, hi=0.85):
    """n colours of monotonically increasing lightness, for an ordered parameter."""
    return [mpl.colormaps[cmap](v) for v in np.linspace(lo, hi, n)]


def panel_label(ax, letter, text='', x=0.0, y=1.03, **kw):
    """Write '(a) text' above the top-left corner of an axes."""
    s = rf'$\bf{{({letter})}}$' + (f' {text}' if text else '')
    return ax.text(x, y, s, transform=ax.transAxes, ha='left', va='bottom',
                   fontsize=kw.pop('fontsize', 9), **kw)


def wigner_image(ax, W, xvec, wlim=W_MAX, cmap='RdBu_r', pvec=None):
    """Draw W(x, p) as a rasterised image with a symmetric colour scale.

    A continuous image (rather than filled contours) avoids the faint banding
    that round-off produces around W = 0 and keeps the vector PDF small.
    pvec: the p grid when it differs from the x grid (a window not centred on the origin).
    """
    pvec = xvec if pvec is None else pvec
    dx, dp = xvec[1] - xvec[0], pvec[1] - pvec[0]
    ext = [xvec[0] - dx / 2, xvec[-1] + dx / 2, pvec[0] - dp / 2, pvec[-1] + dp / 2]
    im = ax.imshow(W, origin='lower', extent=ext, cmap=cmap, vmin=-wlim, vmax=wlim,
                   interpolation='bilinear', rasterized=True)
    ax.set_aspect('equal')
    return im


def time_guides(ax, t_r, labels=False, y=1.02):
    """Vertical guides at t_r/2 (dotted) and t_r (dashed); optional labels on top."""
    ax.axvline(0.5 * t_r, **GUIDE_HALF)
    ax.axvline(t_r, **GUIDE_FULL)
    if labels:
        tr = ax.get_xaxis_transform()
        ax.text(0.5 * t_r, y, r'$t_r/2$', transform=tr, ha='center', va='bottom', fontsize=8,
                color=GUIDE_HALF['color'])
        ax.text(t_r, y, r'$t_r$', transform=tr, ha='center', va='bottom', fontsize=8,
                color=GUIDE_FULL['color'])


def save(fig, stem, png_dpi=300, close=True, copy_to=()):
    """Save <stem>.png and <stem>.pdf at the designed figure size.

    copy_to: extra directories that keep a copy of the same figure (the
    manuscript's supplementary folder mirrors four repository figures).
    """
    import shutil
    stem = Path(stem)
    fig.savefig(stem.with_suffix('.png'), dpi=png_dpi)
    # no creation date in the PDF, so re-running a script reproduces the file
    fig.savefig(stem.with_suffix('.pdf'), dpi=png_dpi, metadata={'CreationDate': None})
    for d in copy_to:
        for ext in ('.png', '.pdf'):
            if (Path(d) / (stem.name + ext)).exists() or ext == '.pdf':
                shutil.copyfile(stem.with_suffix(ext), Path(d) / (stem.name + ext))
    if close:
        plt.close(fig)
    print(f"Saved: {stem.name} (.png, .pdf)")


# ---------------------------------------------------------------------------
# Animation helpers
# ---------------------------------------------------------------------------
def fig_to_frame(fig, dpi=100):
    """Render a figure to a PIL image and close it."""
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, facecolor='white')
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert('RGB')


def save_gif(frames, path, durations_ms):
    """Write a looping GIF.  Durations are in MILLISECONDS per frame.

    (imageio >= 2.28 also expects milliseconds; passing seconds there writes
    0 ms delays and silently drops every pause, which is why this goes through
    Pillow with the unit made explicit.)
    """
    if np.isscalar(durations_ms):
        durations_ms = [durations_ms] * len(frames)
    assert len(durations_ms) == len(frames)
    assert min(durations_ms) >= 20, "GIF frame delays are in milliseconds"
    frames[0].save(path, save_all=True, append_images=frames[1:],
                   duration=[int(d) for d in durations_ms], loop=0, disposal=1)
    print(f"Saved: {Path(path).name}  ({len(frames)} frames, {sum(durations_ms) / 1000:.1f} s)")
