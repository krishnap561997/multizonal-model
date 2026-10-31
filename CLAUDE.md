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
| `room5x5/` | fluid-only spin-up (baseline, 40° vanes) | `userParam05 = 0`, no `startFrom`, `endTime = 600` |
| `room5x5_v2/` | fluid-only spin-up with 60° vanes (jets aimed more at the floor for the smaller room) | `userParam07 = 60`, `endTime = 100`; inlet ≈ 1.55 m/s |
| `room5x5_particles/` | particle stage, restarts from `room0.ini` = last field of `room5x5` | `userParam05 = 1`, 5 M particles, re-seeding every 1 s |

`base-case/`, `base-case2/` are the earlier failed attempts; `docs/base-case-to-room5x5.md` explains
why they failed and what room5x5 changed. Each `room5x5*` folder also has `test.batch` (hpg-dev,
32 ranks, 30 min) for 200-step smoke tests.

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
export NEK_SOURCE_ROOT=$HOME/codes/Nek5000   # Nek5000 master d083ef3 (the scripts default to this)
```
- Source clones live in the **home directory**, never in /blue: Nek5000 at `~/codes/Nek5000`,
  each case's ppiclF at `~/codes/ppiclF_<case folder name>`. The case scripts default to these paths.
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
- Quick tests: `sbatch test.batch` (hpg-dev, `--qos=bala1s`, 32 ranks). Before it, copy `room.par`
  to `room.par.prod` and set `stopAt = numSteps`, `numSteps = 200`; restore afterwards.
  hpg-dev may not accept the burst QOS. A 200-step fluid test takes about 1.5 min (0.41 s/step on 32 ranks).
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
9. **Submitting from an interactive job** (e.g. the VS Code tunnel) leaks its `SLURM_*` step
   variables into the batch job, even with `--export=NONE`. srun then dies in seconds with
   `SLURM_MEM_PER_CPU ... SLURM_MEM_PER_NODE are mutually exclusive` or `Unable to satisfy cpu
   bind request`. Keep the `unset ${!SLURM_CPU_BIND*} ...` block before `srun` in every batch script.

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

## Post-processing
- `tools/nek_slices.py <file with coords> <field file> <outdir>`: |u| contours on the planes
  y = 1.6 m (horizontal, mid-height) and x = 2.5 m (vertical, through the outlet). Only the first
  field file of a run (`room0.f00001`) holds coordinates, so pass it as the coord file.
- System python3 has no matplotlib. A working env (numpy, scipy, matplotlib) is
  `/blue/bala1s/krishnap.kalivel/TLFHydrodynamics/thinfilm-closure-ml/software/envs/thinfilm-ml-py311/bin/python`
  (another project's env: use it, don't modify it).
- PNG output is gitignored; plots go in `<case>/plots/`.

## Current status (update this section as work progresses)
_Last updated 2026-10-31._
- **room5x5** (40°): hpg-dev 200-step test passed at production size (lx1 = 6, 32 ranks), job 48427833.
  Max CFL 0.142, no NaNs, 13.3 GB total memory. No long run yet.
- **room5x5_v2** (60°): built; the bala1s-b run to t = 100 s (job 48429339) was cancelled by Krishna
  because `room.usr` has **no time-averaging**. No output exists.
- **room5x5_particles**: hpg-dev 200-step test passed, job 48433337 (2.5 M particles, 32 ranks,
  restarted from the t = 0.16 s room5x5 test field). Re-seeding worked (1,339 and 1,430 particles
  replaced at t = 0.08 and 0.16 s), max CFL 0.143, no NaNs, 1.39 s/step, 22.2 GB total memory;
  each par file is about 290 MB at 2.5 M particles. For the test, userParam01 = 100 and
  userParam09 = 0.08; production `room.par` is restored. Test output (`par0000*.vtu`, `parout*`,
  `logfile_*`) is still in the folder. A real particle run needs a spun-up `room0.ini`, which doesn't exist yet.
- **Time-averaging is missing** from all room5x5 cases (and from the 2022 v30/v36 cases). `room.usr`
  only declares unused `uavg/vavg/wavg` arrays. `base-case/room.usr` (userchk, ~lines 300–560) has a full
  `avg1/avg2/avg3` + `outpost2` block (files `av1`..`av6`) but averages from step 1 with no reset.
  Proposed next step, awaiting Krishna's go-ahead: new folder `room5x5_v3` (from v2) with running means of
  u, v, w and u², v², w², uv, uw, vw, an averaging start time from `room.par`, output at the field write
  interval. Open question: average within one long run after ~100 s, or spin up 100 s then restart and average.
- Plots so far: `room5x5/plots/slice_*.png` from the 200-step test (t = 0.16 s, flow not developed).

## Conventions
- Keep the 2022 reference cases (`room_one_v30`, `room_one_v36`, Intel 2018 era) read-only.
  They are the physics reference, not build templates.
- Before a production run, make every change in a new, versioned case folder,
  e.g. `room5x5_v2/`.
- Record decisions and results in the case README so Krishna and Sami can follow them.
- Open decisions are for Krishna: the ACH choice for the 5×5 room (inlet Re is 4× lower at the same
  ACH), turning inlet injection back on (v36 style), and periodic boundaries for the multi-room cases.
