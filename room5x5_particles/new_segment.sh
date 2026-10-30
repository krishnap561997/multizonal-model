#!/bin/bash
# Set up the next restart segment of a particle run in its own folder.
#
#   ./new_segment.sh <segNN> <room0.f000XX> <parYYYYY.vtu>
#
# Run it in the folder of the segment that just finished. Pick a fluid field
# file and a par file written at the SAME time: with the defaults
# (writeInterval = 6250 = 5 x IOPART) every room0.f file has a matching par
# file. Check the TIME in the par header:  head -c 300 parYYYYY.vtu
#
# Each segment gets its own folder because ppiclF numbers its output files
# (par/parout/bin) from 00001 in every new job and would overwrite the old
# ones. Note: the first par file of a new segment (par00001.vtu) is a copy
# of parRST.vtu with TIME = 0 in its header - skip it in post-processing.
set -e
SEG=$1; FLD=$2; PAR=$3
if [ -z "$SEG" ] || [ ! -f "$FLD" ] || [ ! -f "$PAR" ]; then
  echo "usage: $0 <segNN> <room0.f000XX> <parYYYYY.vtu>"; exit 1
fi
SRC=$(pwd)
mkdir "$SEG"
cp "$FLD" "$SEG/room0.rst"
cp "$PAR" "$SEG/parRST.vtu"
for f in nek5000 room.re2 room.ma2 room.batch; do cp "$SRC/$f" "$SEG/"; done
sed -e 's/^startFrom *=.*/startFrom = room0.rst/' \
    -e 's/^userParam10 *= *[0-9.]*/userParam10 = 1  /' "$SRC/room.par" > "$SEG/room.par"
echo "Segment $SEG ready:"
grep -E "^startFrom|^userParam10" "$SEG/room.par"
echo "cd $SEG && sbatch room.batch"
