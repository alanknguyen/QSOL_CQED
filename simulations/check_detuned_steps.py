# check_detuned_steps.py
# Step-by-step checks for paper/notes/detuned_derivation.tex: each intermediate
# statement of the derivation (the two-branch form with its phases, the collapse
# envelope, the covariance matrix of a sheared component, the size of the cubic
# term) is compared with the exact amplitudes of the detuned Jaynes-Cummings model.
# Prints only; runs in a few seconds.  The end results are tested in sim_detuned_shear.py.
#
# Nguyen Khoi Nguyen (Alan), Boston University
import numpy as np
from scipy.stats import poisson
g = 1.0; PI = np.pi
def grid(nb):
    n = np.arange(int(nb + 14*np.sqrt(nb) + 40)); return n, np.sqrt(poisson.pmf(n, nb))
def exact(nb, x, phi):
    n, c = grid(nb); Delta = 2*g*np.sqrt(nb)*x; Ob = np.sqrt(Delta**2 + 4*g*g*nb); t = phi*Ob/g**2
    O = np.sqrt(Delta**2 + 4*g*g*(n+1))
    ce = c*(np.cos(O*t/2) - 1j*(Delta/O)*np.sin(O*t/2))
    cg = np.zeros(len(n), complex); cg[1:] = (-1j*c*(2*g*np.sqrt(n+1)/O)*np.sin(O*t/2))[:-1]
    return n, c, O, t, Ob, ce, cg
print("(a) overlap of the exact state with  cos(th/2)|Phi_->(cos(th/2)|e> + e^{i phi} sin(th/2)|g>) + sin(th/2)|Phi_+>(sin(th/2)|e> - e^{-i phi} cos(th/2)|g>)")
print("    with exact branch states |Phi_+-> = sum c_n exp(+-i Omega_n t/2)|n>;  1 - |overlap|^2:")
for nb in (25, 100, 400):
    row = []
    for x in (0.0, 0.5, 1.0, 2.0):
        for phi in (PI/4, PI/2):
            n, c, O, t, Ob, ce, cg = exact(nb, x, phi)
            th = np.arctan2(2*g*np.sqrt(nb), 2*g*np.sqrt(nb)*x)       # mixing angle at n_bar
            Pm, Pp = c*np.exp(-0.5j*O*t), c*np.exp(0.5j*O*t)
            me = np.cos(th/2)**2*Pm + np.sin(th/2)**2*Pp
            mg = np.cos(th/2)*np.sin(th/2)*(np.exp(1j*phi)*Pm - np.exp(-1j*phi)*Pp)
            ov = np.vdot(me, ce) + np.vdot(mg, cg)
            wrong = np.vdot(me, ce) + np.vdot(np.cos(th/2)*np.sin(th/2)*(np.exp(-1j*phi)*Pm - np.exp(1j*phi)*Pp), cg)
            row.append(f"{1-abs(ov)**2:.1e}({1-abs(wrong)**2:.2f})")
    print(f"    n={nb:3d}: " + " ".join(row) + "   [in brackets: same with the sign of phi flipped]")
print("(b) collapse: |sum p_n exp(i Omega_n t)| against exp(-g^2 t^2 sin^2(th)/2), at g t sin(th) = 0.5, 1, 1.5, 2")
for nb in (100, 400):
    for x in (0.0, 1.0, 2.0):
        n, c = grid(nb); Delta = 2*g*np.sqrt(nb)*x; sth = 1/np.sqrt(1+x*x); O = np.sqrt(Delta**2 + 4*g*g*(n+1))
        out = [f"{abs(np.sum(c*c*np.exp(1j*O*(u/sth)))):.4f}/{np.exp(-u*u/2):.4f}" for u in (0.5, 1.0, 1.5, 2.0)]
        print(f"    n={nb:3d} x={x:3.1f}: " + "  ".join(out))
print("(c) covariance of |Phi_+> rotated back by -phi, in (x, p) = (radial, tangential): [Vxx, Vxp, Vpp] against [1/2, -s/2, (1+s^2)/2]; centre angle/phi")
for nb in (100, 1600):
    for x in (0.0, 1.0):
        for phi in (PI/4, PI/2):
            n, c, O, t, Ob, ce, cg = exact(nb, x, phi)
            psi = c*np.exp(0.5j*O*t)
            a1 = np.sum(np.conj(psi[:-1])*psi[1:]*np.sqrt(n[1:])); ang = np.angle(a1)
            psi = psi*np.exp(-1j*ang*n)                              # rotate the centre back onto the x axis
            a1 = np.sum(np.conj(psi[:-1])*psi[1:]*np.sqrt(n[1:])); a2 = np.sum(np.conj(psi[:-2])*psi[2:]*np.sqrt(n[2:]*(n[2:]-1))); nn = np.sum(n*abs(psi)**2)
            da2 = a2 - a1*a1; dn = nn - abs(a1)**2
            Vxx = 0.5 + dn + da2.real; Vpp = 0.5 + dn - da2.real; Vxp = da2.imag
            s = phi/(1+x*x)
            print(f"    n={nb:4d} x={x:.0f} phi={phi/PI:.2f}pi: [{Vxx:.3f}, {Vxp:+.3f}, {Vpp:.3f}] vs [0.500, {-s/2:+.3f}, {(1+s*s)/2:.3f}]   angle/phi = {ang/phi:+.4f}")
print("(d) cubic term: phase phi sin^4(th)/(8 sqrt n) at t_r/2 on resonance:", [f"n={nb}: {PI/2/(8*np.sqrt(nb)):.3f}" for nb in (10, 100, 1000)])
