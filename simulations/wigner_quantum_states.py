# wigner_quantum_states.py
# Static Wigner functions W(x,p) of textbook single-mode states:
# Fock |0>, |1>, |6>, |7>; squeezed vacuum (r = 1, phi = 0 and pi/2);
# coherent |alpha = 2> and |alpha = 2 + 2i>.
#
# QuTiP convention: a = (x + i p)/sqrt(2), so |alpha> is centred at
# (x, p) = sqrt(2) (Re alpha, Im alpha), int W dx dp = 1 and |W| <= 1/pi.
#
# Generates (paper/figures/): wigner_functions_2d, wigner_functions_3d
#
# Nguyen Khoi Nguyen (Alan), Boston University

import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, coherent, squeeze, wigner, ket2dm

from plotstyle import use_style, save, panel_label, wigner_image, FULL_W, W_MAX, PAPER_FIG
from utils import wigner_negativity

use_style()

N = 60
xvec = np.linspace(-6, 6, 301)
X, P = np.meshgrid(xvec, xvec)

states = [
    (r'$|0\rangle$',                       basis(N, 0)),
    (r'$|1\rangle$',                       basis(N, 1)),
    (r'$|6\rangle$',                       basis(N, 6)),
    (r'$|7\rangle$',                       basis(N, 7)),
    (r'sq. vac., $\phi{=}0$',              squeeze(N, 1.0) * basis(N, 0)),
    (r'sq. vac., $\phi{=}\pi/2$',          squeeze(N, 1.0 * np.exp(1j * np.pi / 2)) * basis(N, 0)),
    (r'$|\alpha = 2\rangle$',              coherent(N, 2.0)),
    (r'$|\alpha = 2{+}2i\rangle$',         coherent(N, 2.0 + 2.0j)),
]
W_all = [(lab, wigner(ket2dm(psi), xvec, xvec)) for lab, psi in states]
for lab, W in W_all:
    print(f"{lab:32s} W in [{W.min():+.3f}, {W.max():+.3f}], delta = {wigner_negativity(W, xvec):.3f}")

# -- 2D maps: one colour scale for all panels, zero contour where W changes sign --
fig, axes = plt.subplots(2, 4, figsize=(FULL_W, 3.3), sharex=True, sharey=True)
for ax, letter, (lab, W) in zip(axes.flat, 'abcdefgh', W_all):
    im = wigner_image(ax, W, xvec)
    if W.min() < -1e-3:
        ax.contour(xvec, xvec, W, levels=[0], colors='k', linewidths=0.4, linestyles='--')
    panel_label(ax, letter, lab, fontsize=8)
    ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4])
for ax in axes[-1]:
    ax.set_xlabel(r'$x$')
for ax in axes[:, 0]:
    ax.set_ylabel(r'$p$')
cb = fig.colorbar(im, ax=axes, shrink=0.85, aspect=28, pad=0.02, ticks=[-0.3, -0.15, 0, 0.15, 0.3])
cb.set_label(r'$W(x,p)$')
save(fig, PAPER_FIG / 'wigner_functions_2d')

# -- 3D surfaces with the same colour scale --
fig = plt.figure(figsize=(FULL_W, 3.5))
for k, (letter, (lab, W)) in enumerate(zip('abcdefgh', W_all)):
    ax = fig.add_subplot(2, 4, k + 1, projection='3d')
    ax.plot_surface(X, P, W, cmap='RdBu_r', vmin=-W_MAX, vmax=W_MAX, linewidth=0,
                    antialiased=True, rstride=2, cstride=2, rasterized=True)
    ax.set_zlim(-W_MAX, W_MAX)
    ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4]); ax.set_zticks([-0.3, 0, 0.3])
    ax.tick_params(labelsize=6, pad=-3)
    ax.set_xlabel(r'$x$', labelpad=-9, fontsize=8); ax.set_ylabel(r'$p$', labelpad=-9, fontsize=8)
    ax.view_init(elev=28, azim=-55)
    ax.set_box_aspect((1, 1, 0.75), zoom=0.86)
    ax.set_title(rf'$\bf{{({letter})}}$ {lab}', fontsize=8, pad=-6)
save(fig, PAPER_FIG / 'wigner_functions_3d')
