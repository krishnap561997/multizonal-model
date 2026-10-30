#!/bin/bash
# Regenerate room.re2 and room.ma2 from room.box.  Re-run after ANY edit of room.box.
set -e
NEK_SOURCE_ROOT=${NEK_SOURCE_ROOT:-/blue/bala1s/$USER/codes/Nek5000}
echo room.box | $NEK_SOURCE_ROOT/bin/genbox > genbox.log
mv box.re2 room.re2
printf "room\n0.01\n" | $NEK_SOURCE_ROOT/bin/genmap > genmap.log
# Conforming mesh check: 30x19x30 room + 4x3x4 outlet sharing its ceiling
# nodes has 19295 unique vertices.  More than that = room and outlet boxes
# do not line up (this was the bug in base-case2/room.box).
grep "done locglob_lexico" genmap.log
head -c 132 room.ma2 | strings
