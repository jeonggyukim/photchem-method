"""Hydrogen density cube of the Wolfire, Hollenbach & McKee (2010) cloud for
athinput.cloud_w10: n(r) = N/(2 sqrt(r^2 + r_c^2)), a mean column N = M(r)/(mu_H pi
r^2) equal at every radius (their eq. 2) softened inside the core radius r_c, cut at
r_out and surrounded by n_amb. Written as raw little-endian float32 indexed [k, j, i],
as <problem>/density_file reads it.

usage: python make_density_w10.py OUT NX [--box L] [--column N] [--rc R] [--rout R]
       [--namb N] [--centre X]
       defaults: 40 pc box, N = 1.5e22 cm^-2, r_c = 1 pc, r_out = 17.5 pc,
       n_amb = 0.05 cm^-3, centre on the cell centre nearest the origin.
"""
import argparse
import numpy as np

PC = 3.0857e18
MH, MU_H, MSUN = 1.6735575e-24, 1.4, 1.989e33


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('nx', type=int)
    ap.add_argument('--box', type=float, default=40.0)       # [pc]
    ap.add_argument('--column', type=float, default=1.5e22)  # [cm^-2]
    ap.add_argument('--rc', type=float, default=1.0)         # [pc]
    ap.add_argument('--rout', type=float, default=17.5)      # [pc]
    ap.add_argument('--namb', type=float, default=0.05)      # [cm^-3]
    ap.add_argument('--centre', type=float, default=None)    # [pc], each axis
    a = ap.parse_args()
    dx = a.box/a.nx
    x = -a.box/2 + (np.arange(a.nx) + 0.5)*dx
    xc = a.centre if a.centre is not None else x[a.nx//2]
    z, y, xx = np.meshgrid(x - xc, x - xc, x - xc, indexing='ij')
    r = np.sqrt(xx**2 + y**2 + z**2)
    n = a.column/(2.0*np.sqrt(r**2 + a.rc**2)*PC)
    n = np.where(r <= a.rout, n, a.namb).astype('<f4')
    n.tofile(a.out)
    mass = n[r <= a.rout].astype(np.float64).sum()*MU_H*MH*(dx*PC)**3/MSUN
    print('%s: %d^3, dx = %.4f pc, centre %.5f pc, n_max = %.4g cm^-3, n(r_out) = %.4g, '
          'cloud mass %.4g Msun' % (a.out, a.nx, dx, xc, n.max(),
                                    a.column/(2*np.hypot(a.rout, a.rc)*PC), mass))


if __name__ == '__main__':
    main()
