# Sourced by every build and run script in tests/ (bash or zsh):
#     source "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)/tests/env.sh"
# It reads tests/machine.sh and defines the functions that build a code, run it and
# record what was run. A build writes build.log and a run writes run.log beside its
# output; photchem_method.provenance reads them when a figure script rewrites its
# reduced file, so every reduced file names the code commit it came from.

PHOTCHEM_ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
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
  commit=$(git -C "$dir" rev-parse HEAD) || return 1
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
