"""A clumpy, turbulent-looking cloud for the 3D ray panels: a lognormal density field
(a Gaussian random field with a power-law spectrum, exponentiated) inside a smooth
envelope, with the optical depth a ray accumulates through it.  Used by
adaptive_rays_3d.py and diffuse_rays_3d.py.
"""
import numpy as np


class Cloud:
    """Density on an m^3 grid filling a cube of side `size`, centred at `centre` with
    envelope radius `radius`; rho is normalized to a mean of one inside the envelope."""

    def __init__(self, size, centre, radius, m=40, slope=-3.0, sigma=1.4, seed=3):
        rng = np.random.default_rng(seed)
        k = np.fft.fftfreq(m)*m
        kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
        kk = np.sqrt(kx**2 + ky**2 + kz**2)
        kk[0, 0, 0] = 1.0
        amp = kk**(0.5*slope)
        amp[0, 0, 0] = 0.0
        phase = np.fft.fftn(rng.standard_normal((m, m, m)))
        g = np.real(np.fft.ifftn(amp*phase))
        g = sigma*(g - g.mean())/g.std()
        self.size, self.m = size, m
        self.dx = size/m
        x = (np.arange(m) + 0.5)*self.dx
        X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
        self.centre = np.asarray(centre, dtype=float)
        r2 = ((X - self.centre[0])**2 + (Y - self.centre[1])**2
              + (Z - self.centre[2])**2)/radius**2
        env = np.exp(-r2**2)
        rho = np.exp(g)*env
        self.rho = rho/rho[r2 < 1].mean()
        self.xyz = (X, Y, Z)

    def at(self, pts):
        """Nearest-cell density at points of shape (N, 3); zero outside the cube."""
        idx = np.floor(np.asarray(pts)/self.dx).astype(int)
        ok = np.all((idx >= 0) & (idx < self.m), axis=1)
        out = np.zeros(len(idx))
        i = idx[ok]
        out[ok] = self.rho[i[:, 0], i[:, 1], i[:, 2]]
        return out

    def tau(self, p, q, kappa, nseg=60):
        """Points along p->q and the optical depth accumulated at each, for opacity
        kappa per unit length at the mean density."""
        s = np.linspace(0.0, 1.0, nseg + 1)
        pts = p[None, :] + s[:, None]*(q - p)[None, :]
        rho = self.at(pts)
        ds = np.linalg.norm(q - p)/nseg
        tau = np.concatenate(([0.0], np.cumsum(kappa*ds*0.5*(rho[1:] + rho[:-1]))))
        return pts, tau

    def draw(self, ax, npts=7000, color='0.25', alpha=0.07, ms=3.0, seed=5):
        """Scatter points with probability proportional to density, jittered within
        their cells, so dense clumps read as dark and the envelope as haze."""
        rng = np.random.default_rng(seed)
        w = self.rho.ravel()
        w = w**1.5   # weight the clumps above the haze
        cells = rng.choice(w.size, size=npts, p=w/w.sum())
        X, Y, Z = (a.ravel()[cells] for a in self.xyz)
        jit = (rng.random((3, npts)) - 0.5)*self.dx
        ax.scatter(X + jit[0], Y + jit[1], Z + jit[2], s=ms, color=color, alpha=alpha,
                   edgecolors='none', depthshade=False)
