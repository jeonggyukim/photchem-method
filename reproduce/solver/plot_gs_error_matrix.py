"""Block Gauss-Seidel error matrix of one substep (Appendix E), CO-formation case of F05.
jacobian_dump.txt -> ../../figures/gs_error_matrix.{pdf,png}.
A = I - hJ with J the species Jacobian at frozen coefficients; the species are put in
update order and split into groups; G = -(D_B + L_B)^{-1} U_B for Gauss-Seidel and
G = -D_B^{-1}(L_B + U_B) for Jacobi. Panels show the relative-error form
|G_ij| x_j / x_i (fractional error left in i per fractional error read from j), which
has the same eigenvalues as G and does not depend on the units of the abundances, for
three schemes at one molecular state, with the spectral radius in the title; the bottom
row gives the spectral radius of each scheme along the run."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, '..', '..', 'figures')
path = os.path.join(HERE, 'jacobian_dump.txt')
with open(path) as f:
    names = f.readline().split(':')[1].split()
    myr = float(f.readline().split()[-1])
d = np.loadtxt(path)
nk, ns = int(d[:, 0].max()) + 1, len(names)
J = np.zeros((nk, ns, ns))
X = np.zeros((nk, ns))
t = np.zeros(nk)
for k, tt, i, j, v in d:
    if j < 0:
        X[int(k), int(i)] = v
    else:
        J[int(k), int(i), int(j)] = v
    t[int(k)] = tt
t_myr = t * myr

ORDER = ['He+', 'Si+', 'H2+', 'H3+', 'H+', 'O+', 'C+', 'CHx', 'OHx', 'CO', 'HCO+', 'H2']
perm = [names.index(s) for s in ORDER]
LABEL = {'He+': 'He$^+$', 'Si+': 'Si$^+$', 'H2+': 'H$_2^+$', 'H3+': 'H$_3^+$', 'H+': 'H$^+$',
         'O+': 'O$^+$', 'C+': 'C$^+$', 'CHx': 'CH$_x$', 'OHx': 'OH$_x$', 'CO': 'CO',
         'HCO+': 'HCO$^+$', 'H2': 'H$_2$'}


def groups_default():
    g, i = [], 0
    while i < ns:
        if ORDER[i] == 'CO':
            g.append([i, i + 1]); i += 2
        else:
            g.append([i]); i += 1
    return g


def error_matrix(A, groups, jacobi=False):
    blk = np.zeros(ns, int)
    for b, members in enumerate(groups):
        blk[members] = b
    same = blk[:, None] == blk[None, :]
    lower = blk[:, None] > blk[None, :]
    upper = blk[:, None] < blk[None, :]
    DB = np.where(same, A, 0.0)
    if jacobi:
        return -np.linalg.solve(DB, np.where(lower | upper, A, 0.0))
    return -np.linalg.solve(DB + np.where(lower, A, 0.0), np.where(upper, A, 0.0))


N = 64
h = 3.0 / N
schemes = [('Gauss-Seidel, CO/HCO$^+$ group (default)', groups_default(), False),
           ('Gauss-Seidel, no group', [[i] for i in range(ns)], False),
           ('Jacobi', [[i] for i in range(ns)], True)]
rho = np.zeros((len(schemes), nk))
Gs = {}
for k in range(nk):
    A = np.eye(ns) - h * J[k][np.ix_(perm, perm)]
    for s, (_, groups, jac) in enumerate(schemes):
        G = error_matrix(A, groups, jac)
        rho[s, k] = np.max(np.abs(np.linalg.eigvals(G)))
        Gs[(s, k)] = G

k_show = nk - 1
fig = plt.figure(figsize=(16, 9.5), layout='constrained')
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1])
norm = LogNorm(1e-4, 1)
for s, (title, _, _) in enumerate(schemes):
    ax = fig.add_subplot(gs[0, s])
    xo = np.maximum(X[k_show][perm], 1e-30)
    G = np.abs(Gs[(s, k_show)]) * xo[None, :] / xo[:, None]
    im = ax.imshow(np.where(G > 1e-12, G, np.nan), norm=norm, cmap='viridis')
    ax.set_xticks(range(ns), [LABEL[x] for x in ORDER], rotation=90, fontsize=9)
    ax.set_yticks(range(ns), [LABEL[x] for x in ORDER], fontsize=9)
    ax.set_xlabel('relative error read at the start ($j$)')
    if s == 0:
        ax.set_ylabel('relative error left after the pass ($i$)')
    ax.set_title(f'{title}\nspectral radius {rho[s, k_show]:.3g}', fontsize=10)
fig.colorbar(im, ax=fig.axes[:3], shrink=0.8, label='$|G_{ij}|\\,x_j/x_i$')
ax = fig.add_subplot(gs[1, :])
for s, (title, _, _) in enumerate(schemes):
    ax.semilogy(t_myr, rho[s], marker='o', ms=4, label=title)
ax.axvline(t_myr[k_show], color='grey', ls=':', lw=1)
ax.set_xlabel('t [Myr]')
ax.set_ylabel('spectral radius of $G$')
ax.grid(alpha=0.25, which='both')
ax.legend(fontsize=9)
fig.suptitle(r'Lag error of one substep, CO formation ($n_{\rm H} = 10^3\,{\rm cm^{-3}}$, '
             rf'cosmic rays only), {N} substeps; matrices at t = {t_myr[k_show]:.2f} Myr',
             fontsize=12)
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(FIG, 'gs_error_matrix.' + ext), dpi=200)
for s, (title, _, _) in enumerate(schemes):
    print(f'{title}: rho at t={t_myr[k_show]:.2f} Myr {rho[s, k_show]:.3g}; max over t '
          f'{rho[s].max():.3g} at {t_myr[rho[s].argmax()]:.2f} Myr')
G = Gs[(0, k_show)]
w, v = np.linalg.eig(G)
top = np.argmax(np.abs(w))
vec = np.abs(v[:, top]) / np.abs(v[:, top]).max()
print('default dominant eigenvector:', ', '.join(f'{ORDER[i]} {vec[i]:.2f}'
                                                  for i in np.argsort(-vec)[:4]))
