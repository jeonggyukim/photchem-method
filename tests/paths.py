"""Locations of the codes and runs the harnesses read, from environment variables.

A harness puts reproduce/ on sys.path and imports this module. A variable is read
only when a harness asks for it, so each harness needs only the ones it uses.
"""
import importlib
import os
import sys

VARIABLES = {
    'TIGRIS_DIR': 'a Tigris checkout (branch rayt-photchem-updates)',
    'ATHENAPP_PDR_DIR': 'an Athena++ checkout with '
                        'tst/regression/data/chem_pdr_static.vtk (athena-pp-pdr1d)',
    'ATHENAK_DIR': 'an AthenaK checkout with the chemistry module (athenak-chem)',
    'PHOTCHEM_RUNS': 'the directory holding the Tigris runs (T3_onezone, T6_multi_ion, '
                     'M5_hii_dtype_ions, M6_rad_snr, T7_cost)',
    'PHOTCHEM_POSTPROC_RUNS': 'the directory holding the post-processing runs (be_net)',
    'PYATHENA_DIR': 'a pyathena checkout (its data/ directory is read)',
}


def env(name):
    """The directory named by the environment variable name."""
    value = os.environ.get(name)
    if not value:
        raise RuntimeError('set %s to %s' % (name, VARIABLES[name]))
    return os.path.expanduser(value)


def tigris(*parts):
    return os.path.join(env('TIGRIS_DIR'), *parts)


def athenapp_pdr(*parts):
    return os.path.join(env('ATHENAPP_PDR_DIR'), *parts)


def athenak(*parts):
    return os.path.join(env('ATHENAK_DIR'), *parts)


def runs(*parts):
    return os.path.join(env('PHOTCHEM_RUNS'), *parts)


def postproc_runs(*parts):
    return os.path.join(env('PHOTCHEM_POSTPROC_RUNS'), *parts)


def pyathena(*parts):
    return os.path.join(env('PYATHENA_DIR'), *parts)


def athena_read(code='TIGRIS_DIR'):
    """athena_read from vis/python of the checkout the variable code names."""
    sys.path.insert(0, os.path.join(env(code), 'vis', 'python'))
    return importlib.import_module('athena_read')
