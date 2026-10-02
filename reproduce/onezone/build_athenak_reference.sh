#!/usr/bin/env bash
# Compile athenak_reference.cpp with the include flags of AthenaK's build-cpu
# (read from its CMake flags.make) and link its Kokkos Serial core library.
set -euo pipefail
K=${ATHENAK_DIR:?set ATHENAK_DIR to an AthenaK checkout (athenak-chem)}
B=$K/build-cpu
F=$B/src/CMakeFiles/athena.dir/flags.make
read -r -a INC <<< "$(sed -n 's/^CXX_INCLUDES = //p' "$F")"
${AK_CXX:-/opt/homebrew/bin/g++-16} -O2 -std=c++17 -DKOKKOS_DEPENDENCE "${INC[@]}" \
  "$(dirname "$0")/athenak_reference.cpp" "$B/kokkos/core/src/libkokkoscore.a" \
  -o "$(dirname "$0")/athenak_reference"
