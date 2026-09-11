#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

if [[ -n "${BSIM_ROOT:-}" ]]; then
  CP="$BSIM_ROOT/dist/build:$BSIM_ROOT/lib/core.jar:$BSIM_ROOT/lib/vecmath.jar:$BSIM_ROOT/lib/objimport.jar"
else
  CP="../../dist/build:../../lib/core.jar:../../lib/vecmath.jar:../../lib/objimport.jar"
fi

echo "Compiling HybridDish..."
javac -encoding UTF-8 -cp "$CP" -d . bsim/BSimHybridDish.java bsim/VoxelAnalyzer.java

if [[ $# -eq 0 ]]; then
  echo
  echo "Compiled. From this folder run for example:"
  echo "  ./compile_and_run.sh config/driven_seed101.properties"
  exit 0
fi

echo "Running HybridDish.BSimHybridDish $*"
java -cp ".:$CP" HybridDish.BSimHybridDish "$@"
