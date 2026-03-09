# wigner_fock_states.py
# Wigner functions W(x,p) for Fock states |n>, with marginal distributions.
#
# Nguyen Khoi Nguyen (Alan)
# Dept. of Electrical and Computer Engineering, Boston University
# Advised by Prof. Luca Dal Negro, EC 585
#
# Generates: wigner_fock_individual.png, wigner_fock_combined.png
#
# W_n(x,p) = (-1)^n / pi * L_n(2(x^2 + p^2)) * exp(-(x^2 + p^2))
# L_n: Laguerre polynomial of order n.
# n = 0 (vacuum): Gaussian, no negativity.
# n >= 1: ring structure with negative regions indicating non-classicality.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import eval_laguerre
import matplotlib.gridspec as gridspec


def wigner_fock(n, x, p):
    """W_n(x,p) for Fock state |n>."""
    r2 = x**2 + p**2
    return ((-1)**n / np.pi) * eval_laguerre(n, 2 * r2) * np.exp(-r2)


def plot_individual(n_values, grid_size=100, x_range=5, p_range=5):
    """3D surface plots with marginal distributions on side walls."""
    num = len(n_values)
    rows = int(np.ceil(num / 2))
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(rows, 2, figure=fig)

    x = np.linspace(-x_range, x_range, grid_size)
    p = np.linspace(-p_range, p_range, grid_size)

    for i, n in enumerate(n_values):
        # Adaptive range for high n
        if n > 10:
            rng = 2 * np.sqrt(n + 1)
            xl = np.linspace(-rng, rng, grid_size)
            pl = np.linspace(-rng, rng, grid_size)
        else:
            xl, pl = x, p
        X, P = np.meshgrid(xl, pl)
        W = wigner_fock(n, X, P)

        row, col = i // 2, i % 2
        ax = fig.add_subplot(gs[row, col], projection='3d')
        surf = ax.plot_surface(X, P, W, cmap=cm.coolwarm,
                               linewidth=0, antialiased=True, alpha=0.8)

        # Marginals: integrate W over one quadrature
        dx = xl[1] - xl[0]
        dp = pl[1] - pl[0]
        x_marginal = np.trapezoid(W, dx=dp, axis=0)  # integrate over p
        p_marginal = np.trapezoid(W, dx=dx, axis=1)  # integrate over x

        max_z = np.max(np.abs(W))
        zbase = -1.2 * max_z

        # x marginal on the back wall (p = p_range)
        ax.plot(xl, np.full_like(xl, pl[-1]), zbase + x_marginal * 0.4 / max(np.max(np.abs(x_marginal)), 1e-10),
                color='black', linewidth=2)

        # p marginal on the left wall (x = -x_range)
        ax.plot(np.full_like(pl, xl[0]), pl, zbase + p_marginal * 0.4 / max(np.max(np.abs(p_marginal)), 1e-10),
                color='black', linewidth=2)

        ax.set_xlabel('x')
        ax.set_ylabel('p')
        ax.set_zlabel('W(x,p)')
        ax.set_title(f'Fock State |{n}>')
        ax.view_init(elev=30, azim=-45)
        ax.set_zlim(-max_z * 1.2, max_z * 1.2)

    fig.tight_layout()
    fig.suptitle('Wigner Functions for Fock States', fontsize=16, y=0.98)
    return fig


def plot_combined(n_values, grid_size=100, x_range=5, p_range=5):
    """Compact 2x4 grid of 3D Wigner surfaces."""
    fig = plt.figure(figsize=(15, 8))

    x = np.linspace(-x_range, x_range, grid_size)
    p = np.linspace(-p_range, p_range, grid_size)
    X, P = np.meshgrid(x, p)

    for i, n in enumerate(n_values):
        W = wigner_fock(n, X, P)
        ax = fig.add_subplot(2, 4, i + 1, projection='3d')
        ax.plot_surface(X, P, W, cmap=cm.coolwarm,
                        linewidth=0, antialiased=True, alpha=0.7)
        ax.set_xlabel('x')
        ax.set_ylabel('p')
        ax.set_zlabel('W(x,p)')
        ax.set_title(f'n = {n}')
        ax.view_init(elev=30, azim=-45)
        wmax = np.max(np.abs(W))
        ax.set_zlim(-wmax * 1.2, wmax * 1.2)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])

    fig.tight_layout()
    fig.suptitle('Wigner Functions for Fock States (n = 0 to 7)', fontsize=16, y=0.98)
    return fig


if __name__ == "__main__":
    states = [0, 1, 2, 3, 4, 5, 6, 7]

    fig1 = plot_individual(states)
    fig1.savefig("wigner_fock_individual.png", dpi=300, bbox_inches='tight')

    fig2 = plot_combined(states)
    fig2.savefig("wigner_fock_combined.png", dpi=300, bbox_inches='tight')

    plt.show()
