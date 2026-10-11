"""Where a reduced file's numbers came from: the runs and code checkouts it was read from.

reduced.save calls describe() on the sources a figure script names, so the header of
every reduced file records them. A run directory is described by its run.log, which
tests/env.sh writes (run_logged), and a code checkout by its git commit. Paths are
written relative to the environment variable that locates them ($PHOTCHEM_RUNS/...),
so no machine path enters a tracked file.
"""
import os
import subprocess

from photchem_method.paths import VARIABLES


def portable(path):
    """path with the directory of a set environment variable replaced by ${NAME}."""
    path = os.path.abspath(path)
    best = None
    for name in VARIABLES:
        value = os.environ.get(name)
        if not value:
            continue
        root = os.path.abspath(os.path.expanduser(value))
        if (path == root or path.startswith(root + os.sep)) and \
                (best is None or len(root) > len(best[1])):
            best = (name, root)
    if best is None:
        return path
    return '${%s}%s' % (best[0], path[len(best[1]):])


def revision(path):
    """'branch=... commit=... tree=clean|dirty' of the git checkout containing path."""
    d = path if os.path.isdir(path) else os.path.dirname(path)

    def git(*args):
        return subprocess.run(['git', '-C', d] + list(args), capture_output=True,
                              text=True, check=True).stdout.strip()

    try:
        commit = git('rev-parse', 'HEAD')
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 'not a git checkout'
    dirty = git('status', '--porcelain', '--untracked-files=no')
    return 'branch=%s commit=%s tree=%s' % (git('rev-parse', '--abbrev-ref', 'HEAD'),
                                           commit, 'dirty' if dirty else 'clean')


def describe(source):
    """One line for a source: a run directory (its run.log) or a file or directory in a
    code checkout (its commit)."""
    log = os.path.join(source, 'run.log') if os.path.isdir(source) else None
    if log and os.path.exists(log):
        with open(log) as f:
            code = [ln.split('code: ', 1)[1].strip() for ln in f
                    if ln.startswith('build code: ')]
        return 'run %s (%s)' % (portable(source), code[0] if code else 'run.log has no '
                                'build record')
    if os.path.isdir(source) and not os.path.isdir(os.path.join(source, '.git')) and \
            revision(source) == 'not a git checkout':
        return ('run %s (no run.log: run before tests/env.sh recorded runs; code '
                'revision unknown)' % portable(source))
    return 'file %s (%s)' % (portable(source), revision(source))


def header(sources):
    """The provenance lines of a reduced file: this repository's commit, then one line
    per source."""
    here = os.path.dirname(os.path.abspath(__file__))
    return ['photchem-method %s' % revision(here)] + [describe(s) for s in sources]
