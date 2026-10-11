"""Read lambda_<ion>.txt and interpolate log10 Lambda_q bilinearly in
(log10 T, log10 n_e).  Inputs outside the grid are clamped to its edges.

    from lambda_table import LambdaTable
    tab = LambdaTable('o_3')
    tab(8000.0, 100.0)          # Lambda in erg cm^3 s^-1
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def read_table(path):
    header = {}
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                break
            key, _, val = line[1:].partition(':')
            header[key.strip()] = val.strip()
    logT = np.array(header['logT_grid'].split(), float)
    logne = np.array(header['logne_grid'].split(), float)
    loglam = np.loadtxt(path, comments='#')
    assert loglam.shape == (logT.size, logne.size)
    return logT, logne, loglam, header


class LambdaTable:
    def __init__(self, ion, directory=HERE):
        self.ion = ion
        path = os.path.join(directory, f'lambda_{ion}.txt')
        self.logT, self.logne, self.loglam, self.header = read_table(path)

    def log_lambda(self, T, ne):
        x = np.clip(np.log10(np.asarray(T, float)),
                    self.logT[0], self.logT[-1])
        y = np.clip(np.log10(np.asarray(ne, float)),
                    self.logne[0], self.logne[-1])
        x, y = np.broadcast_arrays(x, y)
        i = np.clip(np.searchsorted(self.logT, x) - 1, 0, self.logT.size - 2)
        j = np.clip(np.searchsorted(self.logne, y) - 1, 0,
                    self.logne.size - 2)
        fx = (x - self.logT[i]) / (self.logT[i + 1] - self.logT[i])
        fy = (y - self.logne[j]) / (self.logne[j + 1] - self.logne[j])
        L = self.loglam
        return ((1 - fx) * (1 - fy) * L[i, j] + fx * (1 - fy) * L[i + 1, j]
                + (1 - fx) * fy * L[i, j + 1] + fx * fy * L[i + 1, j + 1])

    def __call__(self, T, ne):
        return 10.0 ** self.log_lambda(T, ne)
