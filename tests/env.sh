# Sourced by every build and run script in tests/ (bash or zsh):
#     source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
# It finds the repository root from its own path, reads tests/machine.sh and defines
#   code_revision DIR             the commit of a checkout and whether its tree is dirty
#   write_build_log OUTDIR DIR    OUTDIR/build.log, plus OUTDIR/code.diff for a dirty tree
#   tigris_build OUTDIR ARGS      build Tigris into OUTDIR/athena
#   athenak_build OUTDIR FLAGS    build AthenaK into OUTDIR/athena
#   run_logged RUNDIR BUILDDIR CMD  run CMD; write run.out, run.err, time.out, run.log
# photchem_method.provenance reads build.log and run.log when a figure script rewrites
# its reduced file, so every reduced file names the code commit it came from.

# The repository root, from this file's own location (bash and zsh name a sourced
# file differently).
if [ -n "${BASH_VERSION:-}" ]; then
  _env_sh=${BASH_SOURCE[0]}
else
  eval '_env_sh=${(%):-%x}'
fi
PHOTCHEM_ROOT=$(cd "$(dirname "$_env_sh")/.." && pwd)
unset _env_sh
if [ ! -f "$PHOTCHEM_ROOT/tests/machine.sh" ]; then
  echo "tests/machine.sh is missing. Copy the preset for this machine, e.g." \
       "cp tests/machines/grammar.sh tests/machine.sh, and edit it if needed" >&2
  exit 1
fi
. "$PHOTCHEM_ROOT/tests/machine.sh"

# code_revision DIR: one line naming the checkout DIR, its commit and whether the
# source differs from that commit: tree=dirty when a tracked file is modified or an
# untracked file exists under src/.
code_revision() {
  local dir=$1 commit branch dirty
  if ! commit=$(git -C "$dir" rev-parse HEAD 2> /dev/null); then
    echo "dir=$dir not a git checkout: revision unknown"
    return 0
  fi
  branch=$(git -C "$dir" rev-parse --abbrev-ref HEAD)
  dirty=clean
  [ -n "$(git -C "$dir" status --porcelain --untracked-files=no)" ] && dirty=dirty
  [ -n "$(git -C "$dir" ls-files --others --exclude-standard -- src)" ] && dirty=dirty
  echo "repo=$(git -C "$dir" config --get remote.origin.url || echo none)" \
       "branch=$branch commit=$commit tree=$dirty"
}

# write_build_log OUTDIR CODE_DIR CONFIGURE_LINE...: build.log in OUTDIR. For a dirty
# tree it also writes OUTDIR/code.diff (git diff HEAD) and lists in build.log the
# untracked files under src/, so the commit plus code.diff plus those files is the
# source that was built.
write_build_log() {
  local out=$1 dir=$2 code
  shift 2
  code=$(code_revision "$dir")
  {
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "code: $code"
    echo "configure: $*"
    echo "compiler: $($CXX --version | head -n 1)"
    echo "machine: $(uname -sm)"
  } > "$out/build.log"
  case $code in
    *tree=dirty*)
      git -C "$dir" diff HEAD > "$out/code.diff"
      echo "diff: code.diff" >> "$out/build.log"
      git -C "$dir" ls-files --others --exclude-standard -- src \
        | sed 's/^/untracked: /' >> "$out/build.log" ;;
  esac
  # Submodules (AthenaK's kokkos) are part of the source that was built.
  if [ -f "$dir/.gitmodules" ]; then
    git -C "$dir" submodule status | sed 's/^/submodule: /' >> "$out/build.log"
  fi
}

# athenak_build OUTDIR CMAKE_FLAGS...: configure $ATHENAK_DIR out of tree in
# OUTDIR/build with CMAKE_FLAGS plus ATHENAK_CMAKE_FLAGS of machine.sh, build it, copy
# the binary to OUTDIR/athena and write OUTDIR/build.log. The checkout is not touched.
athenak_build() {
  local out=$1
  shift
  local K=${ATHENAK_DIR:?set ATHENAK_DIR to an AthenaK checkout with the chemistry module}
  mkdir -p "$out"
  ("$CMAKE" -S "$K" -B "$out/build" -D CMAKE_BUILD_TYPE=Release \
       "${ATHENAK_CMAKE_FLAGS[@]}" "$@" > "$out/configure.out" 2>&1 \
     && "$CMAKE" --build "$out/build" -j "$MAKE_JOBS" > "$out/build.out" 2>&1) \
    || { echo "build failed; see $out/configure.out and $out/build.out" >&2; return 1; }
  /bin/cp -f "$out/build/src/athena" "$out/athena"
  write_build_log "$out" "$K" cmake "${ATHENAK_CMAKE_FLAGS[@]}" "$@"
}

# tigris_build OUTDIR CONFIGURE_ARGS...: configure and build $TIGRIS_DIR with
# CONFIGURE_ARGS plus the library flags of machine.sh, copy bin/athena to OUTDIR and
# write OUTDIR/build.log.
tigris_build() {
  local out=$1
  shift
  local T=${TIGRIS_DIR:?set TIGRIS_DIR to a Tigris checkout}
  mkdir -p "$out"
  (cd "$T" && python3 configure.py "$@" "${TIGRIS_LIB_FLAGS[@]}" > "$out/configure.out" \
     && make clean > /dev/null && make -j"$MAKE_JOBS" > "$out/build.out" 2>&1) \
    || { echo "build failed; see $out/build.out" >&2; return 1; }
  /bin/cp -f "$T/bin/athena" "$out/athena"
  write_build_log "$out" "$T" "$@" "${TIGRIS_LIB_FLAGS[@]}"
}

_now() { "$PYTHON" -c 'import time; print("%.2f" % time.time())'; }

# run_logged RUNDIR BUILDDIR COMMAND...: run COMMAND in RUNDIR (stdout to run.out,
# stderr to run.err) and write RUNDIR/run.log with the build.log of BUILDDIR, the
# command line, the host, the exit status and the wall time. The wall time also goes
# to time.out as "real SECONDS", the format of /usr/bin/time -p that the figure scripts
# read. COMMAND may be a shell function:
#     run_logged "$W/run" "$W" mpi_launch 4 ../athena -i athinput.local
run_logged() {
  local run=$1 build=$2 status start wall
  shift 2
  mkdir -p "$run"
  {
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    sed 's/^/build /' "$build/build.log"
    echo "command: $*"
    echo "host: $(hostname)"
  } > "$run/run.log"
  start=$(_now)
  (cd "$run" && "$@" > run.out 2> run.err)
  status=$?
  wall=$("$PYTHON" -c "print('%.2f' % ($(_now) - $start))")
  echo "real $wall" > "$run/time.out"
  { echo "exit: $status"; echo "wall_s: $wall"; } >> "$run/run.log"
  return $status
}
