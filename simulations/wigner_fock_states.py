# wigner_fock_states.py
# Wigner functions W(x,p) of the Fock states |n>, n = 0 ... 7, from the closed form
#
#   W_n(x,p) = (-1)^n / pi * L_n(2(x^2 + p^2)) * exp(-(x^2 + p^2)),
#
# with L_n the Laguerre polynomial (convention a = (x + ip)/sqrt(2)).
# n = 0 (vacuum) is a Gaussian with no negativity; for n >= 1 the rings
# alternate in sign and W_n(0,0) = (-1)^n / pi.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Began as a course project in EC 585, taught by Prof. Luca Dal Negro
#
# Generates (paper/figures/): wigner_fock_combined.png, wigner_fock_individual.png

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.special import eval_laguerre

from plotstyle import use_style, FULL_W, W_MAX, PAPER_FIG

use_style()


def wigner_fock(n, x, p):
    """W_n(x,p) for the Fock state |n>."""
    r2 = x ** 2 + p ** 2
    return ((-1) ** n / np.pi) * eval_laguerre(n, 2 * r2) * np.exp(-r2)


def surface(ax, X, P, W):
    # symmetric colour scale: W = 0 is the same (white) in every panel
    ax.plot_surface(X, P, W, cmap='RdBu_r', vmin=-W_MAX, vmax=W_MAX, linewidth=0,
                    antialiased=True, rstride=1, cstride=1, rasterized=True)
    ax.set_zlim(-W_MAX, W_MAX)
    ax.view_init(elev=28, azim=-55)


def plot_combined(n_values, grid_size=121, x_range=4.5):
    """Compact 2 x 4 grid of surfaces."""
    x = np.linspace(-x_range, x_range, grid_size)
    X, P = np.meshgrid(x, x)
    fig = plt.figure(figsize=(FULL_W, 3.5))
    for i, n in enumerate(n_values):
        ax = fig.add_subplot(2, 4, i + 1, projection='3d')
        surface(ax, X, P, wigner_fock(n, X, P))
        ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4]); ax.set_zticks([-0.3, 0, 0.3])
        ax.tick_params(labelsize=6, pad=-3)
        ax.set_xlabel(r'$x$', labelpad=-9, fontsize=8); ax.set_ylabel(r'$p$', labelpad=-9, fontsize=8)
        ax.set_box_aspect((1, 1, 0.75), zoom=0.86)
        ax.set_title(rf'$|{n}\rangle$', fontsize=9, pad=-6)
    return fig


def plot_individual(n_values, grid_size=121, x_range=4.5):
    """Surfaces with the two marginal distributions drawn on the side walls."""
    x = np.linspace(-x_range, x_range, grid_size)
    X, P = np.meshgrid(x, x)
    rows = int(np.ceil(len(n_values) / 2))
    fig = plt.figure(figsize=(FULL_W, 2.6 * rows))
    for i, n in enumerate(n_values):
        W = wigner_fock(n, X, P)
        ax = fig.add_subplot(rows, 2, i + 1, projection='3d')
        surface(ax, X, P, W)
        # marginals |psi_n(x)|^2 and |psi_n(p)|^2 (equal for a Fock state), rescaled for display
        marg_x = trapezoid(W, x, axis=0)          # integrate over p
        marg_p = trapezoid(W, x, axis=1)          # integrate over x
        scale = 0.9 * W_MAX / marg_x.max()
        ax.plot(x, np.full_like(x, x[-1]), -W_MAX + scale * marg_x, color='k', lw=0.8)
        ax.plot(np.full_like(x, x[0]), x, -W_MAX + scale * marg_p, color='k', lw=0.8)
        ax.set_xlabel(r'$x$', labelpad=-6); ax.set_ylabel(r'$p$', labelpad=-6)
        ax.set_zticks([-0.3, 0, 0.3])
        ax.tick_params(labelsize=7, pad=-2)
        ax.set_title(rf'$|{n}\rangle$', fontsize=9, pad=-4)
    return fig


if __name__ == "__main__":
    states = [0, 1, 2, 3, 4, 5, 6, 7]
    fig = plot_combined(states)
    fig.savefig(PAPER_FIG / 'wigner_fock_combined.png', dpi=300)
    plt.close(fig)
    fig = plot_individual(states)
    fig.savefig(PAPER_FIG / 'wigner_fock_individual.png', dpi=250)
    plt.close(fig)
    print("Saved: wigner_fock_combined.png, wigner_fock_individual.png")
