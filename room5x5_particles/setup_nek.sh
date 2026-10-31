#!/bin/bash
# One-time: clone Nek5000 (pinned commit) and build genbox/genmap with the
# same toolchain used for the case.  Run on a HiPerGator login node.
set -e
module purge
module load gcc/14.2.0 openmpi/5.0.10

NEK_SOURCE_ROOT=${NEK_SOURCE_ROOT:-$HOME/codes/Nek5000}
NEK_COMMIT=d083ef340f32fd69cc551a07c7fa56a64c743c05   # master, 2026-10-03 (tested)

if [ ! -d "$NEK_SOURCE_ROOT" ]; then
  mkdir -p "$(dirname "$NEK_SOURCE_ROOT")"
  git clone https://github.com/Nek5000/Nek5000 "$NEK_SOURCE_ROOT"
  (cd "$NEK_SOURCE_ROOT" && git checkout -q $NEK_COMMIT)
fi

cd "$NEK_SOURCE_ROOT/tools"
FC=gfortran CC=gcc ./maketools genbox genmap
echo "Nek5000 ready at $NEK_SOURCE_ROOT  (tools in $NEK_SOURCE_ROOT/bin)"
echo "Add to ~/.bashrc:  export NEK_SOURCE_ROOT=$NEK_SOURCE_ROOT ; export PATH=\$NEK_SOURCE_ROOT/bin:\$PATH"
