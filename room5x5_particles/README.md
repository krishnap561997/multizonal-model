# room5x5_particles — particle stage of the 5 × 5 × 3.2 m room

This is the particle stage of the statistically overloaded run in the 5 × 5 room. It uses the same code, mesh and toolchain as `../room5x5`, which is the fluid-only spin-up. It starts from that run's final flow field.

## What's in it (beyond `../room5x5`)

- **Particles on by default.** `userParam05 = 1`, and `room.par` starts from `room0.ini`.
- **Seeding at t = 0.** 5 M particles in 10 sizes (0.2–50 µm diameter) are placed uniformly in the room. That is the same number density as the 20 M particles in the 2022 10×10 room.
- **Re-seeding** (`userParam09`, default every 1 s):
  - Particles that have left the room (outlet, floor, walls, ceiling) are replaced by new ones. The new particles have the same size, are placed uniformly at random, and start at rest.
  - So each size keeps its initial count, which is the overloading approach of Krishnaprasad et al. 2025.
  - Every particle carries a tag (rank, birth cycle, index): `tag01..03` in the par/parout files. Birth time = birth cycle × dt.
  - Set `userParam09 = 0` to switch re-seeding off. That gives the 2022 behaviour, where the initial particles simply decay.
- **Particle restarts** (`userParam10`, plus `new_segment.sh`): a long particle run can be split into jobs.
  - Particles carry over with their velocities, Langevin state and tags.
  - The cycle numbers continue across jobs, so birth times stay consistent.
- Exits are written to `parout*.vtu` as in 2022. The exit location code is in `iprop_out(5)`: 1 floor, 2/4 x-walls, 3 ceiling, 5/6 z-walls, 7 outlet.

Inlet injection (the v36 `my_place_particle` routine) is still off, as in v30. Periodic boundaries in ppiclF are not needed for a single walled room; they come with the multi-room cases.

**Tested** in a cloud workspace with lx1 = 4, 4 ranks and 40k particles:
- seeding and re-seeding: the counts per size stay constant and the tags are unique, with new particles tagged with their birth cycle;
- a restart from `parRST.vtu`: positions identical, tags preserved, cycle numbers continued.

The production settings (lx1 = 6, 64 ranks, 5 M particles) have not been run yet.

## Steps

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.10
export NEK_SOURCE_ROOT=/blue/bala1s/$USER/codes/Nek5000

./build_ppiclf.sh
./make_mesh.sh                               # or copy room.re2/room.ma2 from ../room5x5
./makenek room
cp ../room5x5/room0.f000NN room0.ini          # last field of the fluid spin-up
sbatch room.batch                             # 64 ranks
```

First do a short hpg-dev test, e.g. `numSteps = 200`, 32 ranks and `userParam08 = 2.5e6`.

**Check the log for:**
- `room5x5: particles per rank per size ...`
- `PPICLF: -Global particles ...`
- `room5x5: re-seeded particles N at t = ...`, once per second.

**Restarting a long run.** The burst QOS (`bala1s-b`) has a 4-day limit; for longer jobs use `--qos=bala1s`.
1. Pick a `room0.f000XX` and a `parYYYYY.vtu` written at the same time. With the defaults, `writeInterval` (6250) is 5 × `userParam01` (1250), so every field file has a matching par file.
2. Run `./new_segment.sh seg02 room0.f000XX parYYYYY.vtu`.
3. `cd seg02 && sbatch room.batch`.

Each segment gets its own folder because ppiclF restarts its output numbering at 00001 in every job. In a new segment, `par00001.vtu` is just a copy of the restart file with TIME = 0 in its header; skip it in post-processing.

## Sizes and cost

- **Particles per rank.** `PPICLF_LPART = 160000` per rank. The run stops at start-up if the particles per rank exceed LPART/2.
  - 5 M particles on 64 ranks is 78k per rank, which is fine.
  - On 32 ranks, use `userParam08 ≤ 2.5e6`, or raise `PPICLF_LPART` and re-run `build_ppiclf.sh` and `makenek`. Memory grows by about 1.2 GB per rank for every 100k of LPART.
- **Disk.** A par file with 5 M particles is about 0.6 GB, written every `userParam01` = 1250 steps (1 s). 1000 s of run time is therefore about 0.6 TB. Raise `userParam01` if `/blue` quota is a concern; parout files are small.
- **Restarts need enough ranks.** `ReadParticleVTU` first loads up to LPART particles onto each of the lowest-numbered ranks. A restart therefore needs `np × LPART ≥` the total number of particles.
