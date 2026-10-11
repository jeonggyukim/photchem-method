"""Named arrays in one text file: the reduced data a figure script plots.

The file holds '#' comment lines, then for each array a line '@ name dim1 dim2 ...'
followed by its values, one row of the array per line (one value per line for a
1-D array). The default format, %.17g, reads back to the same float64.
"""
import numpy as np


def save(path, header, arrays, fmts=None):
    fmts = fmts or {}
    with open(path, 'w') as f:
        for line in header.splitlines():
            f.write('# %s\n' % line)
        for name, value in arrays.items():
            a = np.asarray(value, dtype=np.float64)
            f.write('@ %s %s\n' % (name, ' '.join(str(n) for n in a.shape)))
            np.savetxt(f, a.reshape(a.shape[0], -1) if a.ndim > 1 else a.reshape(-1, 1),
                       fmt=fmts.get(name, '%.17g'))


def load(path):
    arrays, name, shape, values = {}, None, None, []

    def close():
        if name is not None:
            arrays[name] = np.array(values, dtype=np.float64).reshape(shape)

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                close()
                tok = line.split()
                name, shape, values = tok[1], tuple(int(n) for n in tok[2:]), []
            else:
                values.extend(float(v) for v in line.split())
    close()
    return arrays
