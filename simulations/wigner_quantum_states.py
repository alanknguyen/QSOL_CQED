# wigner_quantum_states.py
# Static Wigner functions W(x,p) of textbook single-mode states:
# Fock |0>, |1>, |6>, |7>; squeezed vacuum (r = 1, phi = 0 and pi/2);
# coherent |alpha = 2> and |alpha = 2 + 2i>.
#
# Source for the manuscript's Figs. 10 (3D surfaces) and 11 (2D contours).
# QuTiP convention: a = (x + i p)/sqrt(2), so |alpha> is centred at
# (x, p) = sqrt(2) (Re alpha, Im alpha) and int W dx dp = 1.
#
# Generates (in paper/figures/): wigner_functions_3d.png, wigner_functions_2d.png
#
# Nguyen Khoi Nguyen (Alan), Boston University

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from qutip import basis, coherent, squeeze, wigner, ket2dm

PAPER_FIG = Path(__file__).resolve().parents[1] / 'paper' / 'figures'
PAPER_FIG.mkdir(parents=True, exist_ok=True)

N = 60
xvec = np.linspace(-5, 5, 201)
X, P = np.meshgrid(xvec, xvec)

states = [
    (r'Fock $|0\rangle$ (vacuum)',              ket2dm(basis(N, 0))),
    (r'Fock $|1\rangle$',                       ket2dm(basis(N, 1))),
    (r'Fock $|6\rangle$',                       ket2dm(basis(N, 6))),
    (r'Fock $|7\rangle$',                       ket2dm(basis(N, 7))),
    (r'Squeezed vacuum ($r=1$, $\phi=0$)',      ket2dm(squeeze(N, 1.0) * basis(N, 0))),
    (r'Squeezed vacuum ($r=1$, $\phi=\pi/2$)',  ket2dm(squeeze(N, 1.0 * np.exp(1j * np.pi / 2)) * basis(N, 0))),
    (r'Coherent $|\alpha = 2\rangle$',          ket2dm(coherent(N, 2.0))),
    (r'Coherent $|\alpha = 2 + 2i\rangle$',     ket2dm(coherent(N, 2.0 + 2.0j))),
]
W_all = [(lab, wigner(rho, xvec, xvec)) for lab, rho in states]

# -- 3D surfaces -------------------------------------------------------------
fig = plt.figure(figsize=(12, 20))
for k, (lab, W) in enumerate(W_all):
    ax = fig.add_subplot(4, 2, k + 1, projection='3d')
    wmax = np.max(np.abs(W))
    ax.plot_surface(X, P, W, cmap='coolwarm', vmin=-wmax, vmax=wmax,
                    linewidth=0, antialiased=True, rstride=2, cstride=2)
    ax.set_xlabel(r'$x$', labelpad=4); ax.set_ylabel(r'$p$', labelpad=4)
    ax.set_zlabel(r'$W$', labelpad=4)
    ax.set_zlim(-wmax * 1.05, wmax * 1.05)
    ax.view_init(elev=28, azim=-50)
    ax.tick_params(labelsize=7)
    ax.set_title(lab + (rf'   ($W_{{\min}} = {W.min():.2f}$)' if W.min() < -1e-3 else '   (no negativity)'),
                 fontsize=10, pad=2)
fig.suptitle('Wigner functions of Fock, squeezed-vacuum and coherent states', fontsize=14, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(PAPER_FIG / 'wigner_functions_3d.png', dpi=200, bbox_inches='tight')
fig.savefig(PAPER_FIG / 'wigner_functions_3d.pdf', bbox_inches='tight')
plt.close(fig)

# -- 2D contour maps -----------------------------------------------------------
fig, axes = plt.subplots(4, 2, figsize=(10, 18))
for k, (lab, W) in enumerate(W_all):
    ax = axes.flat[k]
    wmax = np.max(np.abs(W))
    im = ax.contourf(xvec, xvec, W, levels=np.linspace(-wmax, wmax, 81), cmap='RdBu_r')
    if W.min() < -1e-3:
        ax.contour(xvec, xvec, W, levels=[0], colors='k', linewidths=0.5, linestyles='--')
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$p$')
    ax.set_title(lab, fontsize=10)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                      ticks=np.linspace(-wmax, wmax, 5), format='%.2f')
    cb.set_label(r'$W(x,p)$', labelpad=6)
fig.suptitle('Wigner functions: negative regions (blue, dashed zero contour) are non-classical',
             fontsize=13, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.985])
fig.savefig(PAPER_FIG / 'wigner_functions_2d.png', dpi=200, bbox_inches='tight')
fig.savefig(PAPER_FIG / 'wigner_functions_2d.pdf', bbox_inches='tight')
plt.close(fig)
print("Saved: wigner_functions_3d, wigner_functions_2d")
