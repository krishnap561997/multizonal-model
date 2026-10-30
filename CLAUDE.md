# CLAUDE.md — P3 Multizonal Aerosol Model (HiPerGator)

Context for Claude Code working on HiPerGator. Put this file in the folder that holds
the case folders, e.g. `/blue/bala1s/<user>/multizonal/` containing `room5x5/` and
`room5x5_particles/`.

## Project in one paragraph
We predict aerosol dispersal across connected indoor zones (rooms, halls, auditoriums)
that each have their own ceiling cassette ventilation unit. A multizonal well-mixed model
(paper in progress) is validated against **statistically overloaded Euler–Lagrange LES**:
Nek5000 (spectral-element LES, dynamic Smagorinsky) solves the airflow. ppiclF tracks
millions of one-way-coupled droplet nuclei (10 sizes, Langevin sub-grid velocity, Stokes
settling) so that source–sink statistics for every pair of locations come from one run
(Salinas et al. 2022, PRF; Krishnaprasad et al. 2025). PI/mentor: Krishna (K. A. Krishnaprasad).
Simulations are run by undergraduate Sami Latouf.

Roadmap: (1) single 5×5×3.2 m room (current); (2) four 5×5 rooms in a row, periodic in x;
(3) 2×2 rooms, periodic in x and z. Later stages need ppiclF periodicity and counters
for particles crossing periodic boundaries / zones.

## Cases
| folder | purpose | key settings |
|---|---|---|
| `room5x5/` | fluid-only spin-up | `userParam05 = 0`, no `startFrom` |
| `room5x5_particles/` | particle stage, restarts from `room0.ini` = last field of `room5x5` | `userParam05 = 1`, 5 M particles, re-seeding every 1 s |

Each folder has a README.md. Read it before changing anything.

Geometry:
- Room box 30×19×30 elements (0 ≤ x,z ≤ 5 m, 0 ≤ y ≤ 3.2 m; **y is vertical**).
- Outlet box 4×3×4 elements (2.2–2.8 m in x and z, y 3.2–3.7 m) on top of the room.
- 17,148 elements in total, polynomial order lx1 = 6.
- Four inlet slots 0.055 × 0.44 m sit around the outlet. The inlet velocity is set in `userbc` from ACH (`userParam06`) and the vane angle (`userParam07`).

## Toolchain (do not mix)
```bash
module purge
module load gcc/14.2.0 openmpi/5.0.10
export NEK_SOURCE_ROOT=/blue/bala1s/$USER/codes/Nek5000   # Nek5000 master d083ef3
```
- Each case builds its **own** ppiclF (`./build_ppiclf.sh`, github.com/dpzwick/ppiclF @ 72483db).
  `ppiclf_user.f` and `PPICLF_USER.h` are compiled *into* `libppiclF.a`. After editing
  them, re-run `./build_ppiclf.sh` and then `./makenek room`.
- Do NOT use the old shared `~/Nek5000` or `~/ppiclF` from 2022, and do not use Intel
  compilers. Mixing compilers caused earlier failures, e.g. gfortran-built
  libppiclF linked into Intel Nek → `_gfortran_*` undefined references.
- `genbox`/`genmap` come from `$NEK_SOURCE_ROOT/bin` (built by `setup_nek.sh`).

## Standard workflow
```bash
./setup_nek.sh        # once per user: clone Nek5000, build genbox/genmap
./build_ppiclf.sh     # per case
./make_mesh.sh        # room.box -> room.re2 + room.ma2; expect 19295 unique vertices
./makenek room        # -> ./nek5000 ; look for "Compilation successful" in build.log
sbatch room.batch     # account bala1s, qos bala1s-b (4-day limit), 1 node, 64 ranks
```
- Quick tests: use `--partition=hpg-dev`, a short `--time`, `numSteps = 200` and 32 ranks.
  hpg-dev may not accept the burst QOS; use `--qos=bala1s` there.
- Never run `nek5000` on a login node. Compiling there is fine.
- Logs: `logfile_<jobid>.out`. Healthy start-up prints the line
  `room5x5: ACH 5.00 vane angle 40.0 deg -> inlet speed 1.7857 m/s`, and the CFL number `C=`
  stays below about 0.6. With particles you also see `PPICLF: -Global particles` and,
  every second, `room5x5: re-seeded particles N`.

## Rules that prevent known failures
1. **No `--distribution=cyclic:cyclic`** in Slurm scripts. The current Nek5000 file reader
   uses one aggregator per node and assumes each node holds a contiguous block of ranks. Cyclic placement gives
   `Error reading .re2 mesh ierr=<64*10>` (error 10 = NEKIO_RDERR_USRBUF).
2. **MPI ranks ≥ `lpmin` (32) in SIZE.** If the element count changes, update `lelg`,
   `lelx/lely/lelz` (room box only) in SIZE, `room.box`, and the geometry block at the top of
   `usrdat2` in `room.usr` together.
3. **Room and outlet meshes must conform.** The room-box x/z node lists must contain the outlet
   box nodes (2.2, 2.35, 2.5, 2.65, 2.8). Check that `make_mesh.sh` reports the expected number of
   unique vertices. If there are more, the boxes don't match.
4. **Regenerate `room.re2`/`room.ma2` after every `room.box` edit.**
5. **Particles per rank**: `PPICLF_LPART = 160000`. `room.usr` aborts if the particles per rank
   exceed LPART/2. 5 M particles needs 64 ranks; on 32 ranks use ≤ 2.5 M. Memory grows about
   1.2 GB/rank per 100k LPART.
6. **Fixed-form Fortran** (`.usr`, `ppiclf_user.f`): code must stay within column 72, and
   continuation lines take a character in column 6. Watch for names that clash with Nek
   common blocks: `xc`, `zc`, `wx`, `wy`, `wz` are already Nek arrays.
7. Particle restarts: run `./new_segment.sh segNN room0.f000XX parYYYYY.vtu` with a field file
   and a par file from the *same* time. Each segment runs in its own folder, because ppiclF restarts
   output numbering at 00001. The first `par00001.vtu` of a segment repeats the restart file
   (TIME = 0) and should be skipped.
8. Disk: each par file with 5 M particles is about 0.6 GB, written every 1 s. Check `/blue` quota
   (`blue_quota`) before long runs.

## Code map (`room.usr`)
- `uservp` / `userf`: eddy viscosity plus a sponge/buffer in the outlet box.
- `userchk`:
  - dynamic Smagorinsky;
  - Langevin wall-damping functions f1–f3 (particles only);
  - inlet random-perturbation array;
  - ppiclF interpolation, integration and `parout` output;
  - re-seeding (`reseed_particles`).
- `userbc`: inlet velocity on the ceiling, using a tanh shape function times the inlet speed derived from ACH.
- `usrdat2`:
  - the geometry block ("USER DEFINED GEOMETRY");
  - room/outlet boundary conditions, with ceiling elements under the outlet found from their coordinates;
  - ppiclF initialisation: fresh seeding, or a restart from `parRST.vtu`.
- `ppiclf_user.f`: particle equations of motion and the Langevin model. Particles are removed at
  walls/outlet, and the exit location code goes into `iprop_out(5)`: 1 floor, 2/4 x-walls,
  3 ceiling, 5/6 z-walls, 7 outlet.
- Particle tags (`tag01..03` in VTU files) = (rank, birth cycle, index). Birth time = cycle × dt.

## Conventions
- Keep the 2022 reference cases (`room_one_v30`, `room_one_v36`, Intel 2018 era) read-only.
  They are the physics reference, not build templates.
- Before a production run, make every change in a new, versioned case folder,
  e.g. `room5x5_v2/`.
- Record decisions and results in the case README so Krishna and Sami can follow them.
- Open decisions are for Krishna: the ACH choice for the 5×5 room (inlet Re is 4× lower at the same
  ACH), turning inlet injection back on (v36 style), and periodic boundaries for the multi-room cases.
