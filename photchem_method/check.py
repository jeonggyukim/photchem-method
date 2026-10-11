"""Check the repository against the conventions in tests/README.md.

    python -m photchem_method.check

Exits 1 and lists each violation as path:line: rule: text. The pre-commit hook in
.githooks/ runs it, so a commit that breaks a rule is refused. The rules, on tracked
files: no figure number in a path or output name; a figure script in tests/ sits in a
fig_<name>/ folder, which holds a plot script and, if that script reads --from-runs, a
reduced file with provenance lines; a figure script (schematics/*.py not starting
with _, tests/*/fig_*/plot*.py) that saves a file also creates the output directory
with os.makedirs or mkdir; no machine path outside tests/machines/; no
sys.path import; every environment variable is listed in paths.VARIABLES; a shell
script that configures a build sources tests/env.sh; every shell script passes
bash -n, or zsh -n when its shebang names zsh.
"""
import os
import re
import shutil
import subprocess
import sys

from photchem_method.paths import VARIABLES

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE = re.compile(r'(\.(py|sh|cpp|hpp|md|in|toml)|/athinput\.[^/]*)$')
MACHINES = 'tests/machines/'
MACHINE_PATH = re.compile(r'/opt/homebrew|/usr/local/Cellar|/Users/\w|/home/\w|'
                          r'~/(Documents|Dropbox|Projects)')
FIGNUM_PATH = re.compile(r'(^|/)F\d\d[a-z]?_')
FIGNUM_OUT = re.compile(r"""['"]F\d\d[a-z]?_\w*['"]""")
SYSPATH = re.compile(r'sys\.path\.insert\(0, str\(.*__file__')
ENV_USE = re.compile(r"""paths\.env\(['"]([A-Z_]+)['"]\)|\$\{([A-Z_]+):\?|"""
                     r"""getenv\("([A-Z_]+)"\)""")
ENV_NAME = re.compile(r'_(DIR|RUNS|EXE)$|^XUVTOP$')
BUILD = re.compile(r'configure\.py|cmake ')
PROVENANCE = '# provenance:'
FIGURE_SCRIPT = re.compile(r'^schematics/[^_/][^/]*\.py$|'
                           r'^tests/[^/]+/fig_[^/]+/plot[^/]*\.py$')
SAVES = re.compile(r'savefig\(|\.save\(')
MAKES_DIR = re.compile(r'os\.makedirs\(|\.mkdir\(')


def tracked():
    out = subprocess.run(['git', 'ls-files'], cwd=ROOT, capture_output=True, text=True,
                         check=True).stdout
    return [p for p in out.splitlines() if os.path.isfile(os.path.join(ROOT, p))]


def lines(path):
    with open(os.path.join(ROOT, path), errors='replace') as f:
        return f.read().splitlines()


def shell_syntax(path, text):
    """'' if `<shell> -n` accepts the script, else the first error line. The shell is
    the one its shebang names (zsh or bash; sh and sourced files without a shebang
    are read by bash). A missing shell is reported, not skipped."""
    shell = 'zsh' if text and 'zsh' in text[0] else 'bash'
    if shutil.which(shell) is None:
        return '%s is not installed, so the syntax was not checked' % shell
    r = subprocess.run([shell, '-n', os.path.join(ROOT, path)], capture_output=True,
                       text=True)
    if r.returncode == 0:
        return ''
    return (r.stderr.strip().splitlines() or ['%s -n failed' % shell])[0]


def check(files):
    found = []

    def bad(path, n, rule, text):
        found.append('%s:%d: %s: %s' % (path, n, rule, text.strip()[:100]))

    fig_dirs = {os.path.dirname(p) for p in files
                if re.search(r'(^|/)fig_[^/]+/', p)}
    for path in files:
        if FIGNUM_PATH.search(path):
            bad(path, 0, 'figure-number-in-path', path)
        if path.startswith('tests/') and re.search(r'(^|/)plot[^/]*\.py$', path) \
                and os.path.dirname(path) not in fig_dirs:
            bad(path, 0, 'figure-script-outside-fig-folder', path)
        if not SOURCE.search(path):
            continue
        text = lines(path)
        if path.endswith('.sh'):
            problem = shell_syntax(path, text)
            if problem:
                bad(path, 0, 'shell-syntax', problem)
        if FIGURE_SCRIPT.search(path) and any(SAVES.search(t) for t in text) \
                and not any(MAKES_DIR.search(t) for t in text):
            bad(path, 0, 'figure-script-without-makedirs', path)
        for n, line in enumerate(text, 1):
            if not path.startswith(MACHINES) and path != 'photchem_method/check.py' \
                    and MACHINE_PATH.search(line):
                bad(path, n, 'machine-path', line)
            if path.endswith('.py') and SYSPATH.search(line):
                bad(path, n, 'sys-path-import', line)
            if path.endswith('.py') and FIGNUM_OUT.search(line) and 'paths.runs' not in line:
                bad(path, n, 'figure-number-in-output', line)
            for m in ENV_USE.finditer(line):
                name = next(g for g in m.groups() if g)
                if ENV_NAME.search(name) and name not in VARIABLES:
                    bad(path, n, 'unlisted-env-variable', name)
        if path.endswith('.sh') and not path.startswith(MACHINES) \
                and path != 'tests/env.sh' and any(BUILD.search(t) for t in text) \
                and not any('tests/env.sh' in t for t in text):
            bad(path, 0, 'build-without-env-sh', path)
    for d in sorted(fig_dirs):
        if not d.startswith('tests/'):
            continue
        here = [p for p in files if os.path.dirname(p) == d]
        plots = [p for p in here if re.search(r'/plot[^/]*\.py$', p)]
        if not plots:
            bad(d, 0, 'fig-folder-without-plot', d)
        reads_runs = any('--from-runs' in t for p in plots for t in lines(p))
        reduced = [p for p in here if re.search(r'/reduced[^/]*\.txt$', p)]
        if reads_runs and not reduced:
            bad(d, 0, 'from-runs-without-reduced-file', d)
        for p in reduced:
            if not any(t.startswith(PROVENANCE) for t in lines(p)[:40]):
                bad(p, 0, 'reduced-file-without-provenance', p)
    return found


def main():
    found = check(tracked())
    for f in found:
        print(f)
    if found:
        print('%d violation(s); see "Conventions" in tests/README.md' % len(found))
        sys.exit(1)


if __name__ == '__main__':
    main()
