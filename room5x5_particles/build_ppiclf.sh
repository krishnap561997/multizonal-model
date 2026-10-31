#!/bin/bash
# Build ppiclF for THIS case.
# ppiclF compiles ppiclf_user.f and PPICLF_USER.h into libppiclF.a, so every
# case needs its own ppiclF build. Re-run this script whenever you change
# ppiclf_user.f or PPICLF_USER.h, then re-run ./makenek room.
#
# Usage:  ./build_ppiclf.sh            (inside the case folder, modules loaded)
set -e
CASE_DIR=$(cd "$(dirname "$0")" && pwd)
PPICLF_ROOT=${PPICLF_ROOT:-$HOME/codes/ppiclF_$(basename "$CASE_DIR")}
PPICLF_COMMIT=72483db            # last commit of github.com/dpzwick/ppiclF (2020-02-08)

if [ ! -d "$PPICLF_ROOT" ]; then
  git clone https://github.com/dpzwick/ppiclF "$PPICLF_ROOT"
  (cd "$PPICLF_ROOT" && git checkout -q $PPICLF_COMMIT)
fi

cp "$CASE_DIR/ppiclf_user.f"  "$PPICLF_ROOT/source/ppiclf_user.f"
cp "$CASE_DIR/PPICLF_USER.h"  "$PPICLF_ROOT/source/PPICLF_USER.h"

cd "$PPICLF_ROOT"
rm -f source/*.o source/libppiclF.a source/ppiclf.f source/PPICLF
# -std=legacy also allows the old argument-type mismatches (gfortran >= 10)
make FC=mpif90 CC=mpicc FFLAGS="-cpp -std=legacy -O2 -mcmodel=medium" \
     2>&1 | tee "$CASE_DIR/build_ppiclf.log"
test -f source/libppiclF.a
echo "ppiclF library: $PPICLF_ROOT/source/libppiclF.a"
