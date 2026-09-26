"""Search for a better species update order on the frozen GOW17 Jacobians of the
CO-formation test (jacobian_dump.txt), with CO and HCO+ kept together as a group.

For an order, at each sample k (substep h = t_end/64), two scores in fractional form:
  rho   spectral radius of G = -(D_B + L_B)^{-1} U_B (persistent lag);
  lag   sum over couplings read old (above the block diagonal) of
        |A_ij| x_j / (A_ii x_i) * delta_j, delta_j the fractional change of x_j in one
        substep (lag carried in one substep).
An order's score is the worst over the samples.  Local search by moving one unit
(a species or the CO/HCO+ group) to another position, from the current solver order
and from random starts."""
import io
import contextlib
import os
import runpy
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(os.path.join(HERE, 'plot_gs_error_matrix.py'))
J, X, names, t_myr, myr = g['J'], g['X'], g['names'], g['t_myr'], g['myr']
ORDER0 = g['ORDER']
N = 64
h = 3.0/N
nk = len(t_myr)

# per sample: A, x, delta in the network's own species order
A_k, x_k, d_k = [], [], []
for k in range(1, nk):
    A_k.append(np.eye(len(names)) - h*J[k])
    x = np.maximum(X[k], 1e-30)
    xp = np.maximum(X[k - 1], 1e-30)
    dt = (t_myr[k] - t_myr[k - 1])/myr
    x_k.append(x)
    d_k.append(np.abs(x - xp)/dt*h/x)


def units_to_perm(units):
    return [names.index(s) for u in units for s in (u if isinstance(u, tuple) else (u,))]


def blocks_of(units):
    out, i = [], 0
    for u in units:
        n = len(u) if isinstance(u, tuple) else 1
        out.append(list(range(i, i + n)))
        i += n
    return out


def score(units):
    perm = units_to_perm(units)
    blocks = blocks_of(units)
    blk = np.zeros(len(perm), int)
    for b, m in enumerate(blocks):
        blk[m] = b
    lower = blk[:, None] > blk[None, :]
    upper = blk[:, None] < blk[None, :]
    same = blk[:, None] == blk[None, :]
    rho_w, lag_w = 0.0, 0.0
    for A, x, d in zip(A_k, x_k, d_k):
        Ap = A[np.ix_(perm, perm)]
        xp = x[perm]
        G = -np.linalg.solve(np.where(same | lower, Ap, 0.0), np.where(upper, Ap, 0.0))
        rho_w = max(rho_w, np.max(np.abs(np.linalg.eigvals(G))))
        C = np.abs(Ap)*xp[None, :]/(np.abs(np.diag(Ap))[:, None]*xp[:, None])
        lag_w = max(lag_w, np.sum(np.where(upper, C*d[perm][None, :], 0.0)))
    return rho_w, lag_w


def units_from_order(order):
    u, i = [], 0
    while i < len(order):
        if order[i] == 'CO':
            u.append(('CO', 'HCO+')); i += 2
        else:
            u.append(order[i]); i += 1
    return u


def label(units):
    return ' '.join('{CO,HCO+}' if isinstance(u, tuple) else u for u in units)


def local_search(units, key):
    best = list(units); best_s = score(best)
    improved = True
    while improved:
        improved = False
        for a in range(len(best)):
            for b in range(len(best)):
                if a == b:
                    continue
                cand = list(best); u = cand.pop(a); cand.insert(b, u)
                s = score(cand)
                if key(s) < key(best_s) - 1e-12:
                    best, best_s, improved = cand, s, True
    return best, best_s


base = units_from_order(ORDER0)
s0 = score(base)
print('current order : rho %.3f  lag %.3f   %s' % (s0[0], s0[1], label(base)))
for name, key in (('rho', lambda s: s[0]), ('lag', lambda s: s[1]),
                  ('rho+lag', lambda s: s[0] + s[1])):
    rng = np.random.default_rng(1)
    starts = [base] + [list(np.array(base, dtype=object)[rng.permutation(len(base))])
                       for _ in range(20)]
    results = [local_search(st, key) for st in starts]
    best, bs = min(results, key=lambda r: key(r[1]))
    print('best by %-7s: rho %.3f  lag %.3f   %s' % (name, bs[0], bs[1], label(best)))
