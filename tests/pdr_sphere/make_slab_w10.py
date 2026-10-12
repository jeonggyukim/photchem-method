"""Hydrogen density of the plane-parallel slab Wolfire, Hollenbach & McKee (2010)
solve, for athinput.slab_w10 (nx1 x 2 x 2, raw little-endian float32 indexed [k, j,
i], as <problem>/density_file reads it).

Their PDR code lights a slab with one beam at 60 degrees to the normal, so a layer at
normal depth z sees the beam through 2 N(z), N(z) the normal column, while its
chemistry uses the density there. A beam along +x reproduces both: put at x the
density f n(z) of the cloud at depth z = f x/2 (the 1/r profile of make_density_w10.py
from the surface r_out inward, mirrored past the centre), so that the column along x
is 2 N(z). f = 1 uses the mean density n_bar; f = 4 is their clump density n_c ~ 4
n_bar (their sec. 5.5), at the same beam column.

usage: python make_slab_w10.py OUT NX LENGTH [--factor F] [--column N] [--rc R]
       [--rout R]   (LENGTH: the box along x [pc]; defaults as make_density_w10.py)
"""
import argparse
import numpy as np

PC = 3.0857e18


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('nx', type=int)
    ap.add_argument('length', type=float)                    # [pc]
    ap.add_argument('--factor', type=float, default=1.0)
    ap.add_argument('--column', type=float, default=1.5e22)  # [cm^-2]
    ap.add_argument('--rc', type=float, default=1.0)         # [pc]
    ap.add_argument('--rout', type=float, default=17.5)      # [pc]
    a = ap.parse_args()
    dx = a.length/a.nx
    x = (np.arange(a.nx) + 0.5)*dx
    z = a.factor*x/2.0
    r = np.abs(a.rout - z)
    n = a.factor*a.column/(2.0*np.sqrt(r**2 + a.rc**2)*PC)
    np.repeat(n.astype('<f4')[None, None, :], 2, axis=1).repeat(2, axis=0).tofile(a.out)
    av_normal = np.sum(n*dx*PC)/2.0/1.9e21    # the beam column is twice the normal
    print('%s: %d cells of %.4f pc, factor %g, n %.4g-%.4g cm^-3, normal A_V through '
          'the box %.2f' % (a.out, a.nx, dx, a.factor, n.min(), n.max(), av_normal))


if __name__ == '__main__':
    main()
