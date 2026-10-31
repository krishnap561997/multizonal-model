# room5x5 — single 5 × 5 × 3.2 m room, statistically overloaded LES

This case uses current Nek5000 and ppiclF to simulate a 5 × 5 × 3.2 m room with a centred 4-way cassette, run on 32 or 64 cores on HiPerGator. It was ported from `room_one_v30` (2022).

Build and run were tested with gcc 14.2, Nek5000 master (`d083ef3`) and ppiclF (`72483db`):

* A fluid-only spin-up ran cleanly.
* A particle restart ran cleanly: 40k particles, no NaNs, and particles leaving through the floor.

These tests used a reduced polynomial order (lx1 = 4) and 4 ranks in a cloud workspace. The production SIZE (lx1 = 6, lpmin = 32) has not been run yet; do that first on hpg-dev.

## Files

| file | purpose |
|---|---|
| `room.box` | genbox input. Room box 30×19×30 (0.169 m, 0.15 m under the outlet) plus outlet box 4×3×4 (0.6 × 0.6 × 0.5 m). The grids conform. |
| `SIZE` | lelx/lely/lelz = 30/19/30, lelg = 17148, lpmin = 32, lx1 = 6 |
| `room.par` | Run controls. userParam05 = particles on/off, 06 = ACH, 07 = vane angle, 08 = total particles |
| `room.usr` | Physics: dynamic Smagorinsky, inlet BC, outlet buffer, Langevin wall functions, ppiclF calls. Geometry is set in `usrdat2` ("USER DEFINED GEOMETRY"). |
| `ppiclf_user.f`, `PPICLF_USER.h` | Particle equations and ppiclF array sizes. These are compiled *into* libppiclF. |
| `setup_nek.sh` | One-time setup: clones Nek5000 and builds genbox/genmap |
| `build_ppiclf.sh` | Builds ppiclF for this case (re-run whenever you edit `ppiclf_user.f` or `PPICLF_USER.h`) |
| `make_mesh.sh` | Makes `room.re2` and `room.ma2` from `room.box` |
| `makenek` | Builds `nek5000` |
| `room.batch` | Slurm script: account bala1s, QOS bala1s-b, 1 node, 64 ranks |

## Steps on HiPerGator

```bash
module purge; module load gcc/14.2.0 openmpi/5.0.10
export NEK_SOURCE_ROOT=/blue/bala1s/$USER/codes/Nek5000

./setup_nek.sh            # once per user
./build_ppiclf.sh         # once per case (and after editing ppiclf files)
./make_mesh.sh            # once per mesh; expect "19295" unique vertices
./makenek room            # builds ./nek5000; check build.log for "Compilation successful"
```

**Short test first:** in `room.par` set `numSteps = 200` (instead of `stopAt`/`endTime`). In `room.batch`, switch to the hpg-dev lines with 32 ranks, then run `sbatch room.batch`. Check `logfile_<jobid>.out` for:

* `room5x5: ACH 5.00 vane angle 40.0 deg -> inlet speed 1.7857 m/s`
* steps advancing, with CFL `C=` below about 0.6.

**Stage 1 — fluid spin-up** (userParam05 = 0, no `startFrom`)
- Run to statistically steady flow, a few hundred seconds of physical time.
- Fields are written to `room0.f000NN`.

**Stage 2 — particles**
1. `cp room0.f000NN room0.ini` (last field file).
2. In `room.par`, uncomment `startFrom = room0.ini time=0` and set `userParam05 = 1`.
3. Set `userParam08`, the total particle count. The default 5 M gives the same number density as the 2022 20 M / 10×10 room.
4. Resubmit.

## Things to know

- **Inlet speed comes from ACH.** `uin = ACH·V/3600 / (4·A_inlet·sin θ)`. This reproduces 7.1428 m/s for the 10×10 room at ACH 5 and 40°. For the 5×5 room at ACH 5 it is 1.79 m/s, which means about 4× lower inlet Re. Change ACH and angle in `room.par` only.
- **Particle memory.** `PPICLF_LPART = 160000` per rank. `room.usr` stops if particles per rank exceed LPART/2. With 5 M particles, 64 ranks is fine (78k/rank); 32 ranks is not (156k/rank). For 32 ranks, either lower `userParam08` or raise `PPICLF_LPART`. Raising LPART grows static memory per rank (about 1 GB per 100k), and you must re-run `build_ppiclf.sh` and `makenek`.
- **Never use `--distribution=cyclic:cyclic`.** The current Nek5000 reader then fails with `Error reading .re2 mesh ierr=…`.
- **Keep one toolchain.** Don't mix compilers in one build. Nek5000, ppiclF and gslib must all be built with the same `module load gcc/14.2.0 openmpi/5.0.10`. Don't use the old shared `~/Nek5000` or `~/ppiclF`.
- If `build_ppiclf.sh` fails at "makeThird", the gslib download was empty. Delete `ppiclF/3rd_party/gslib/v1.0.5.tar.gz` and re-run.

## Changes from room_one_v30 (2022)

1. Geometry moved to the top of `usrdat2`. The outlet ceiling elements are now found from their coordinates; no element indices are hard-coded.
2. Inlet speed is computed from ACH and vane angle. `ux` is now zeroed in `userbc`; before, it kept its value from the previous call.
3. `ran2` and `ppiclf_ran` now SAVE `idum1`. Under gfortran it is otherwise undefined between calls. The Langevin generator is also now seeded per rank; it was never seeded before, so every rank drew the same sequence.
4. Loops over `lelv` changed to `nelv`. Outlet-box elements are skipped in the wall-function arrays, which were previously indexed out of bounds.
5. Changes in `userchk`:
   - `gop` now uses its own work array; before, it overwrote Nek's mesh-velocity arrays.
   - Wall-distance arrays are SAVEd.
   - A missing comma in a declaration was fixed.
6. Particles can be switched on/off with userParam05; fluid-only runs skip ppiclF entirely. Inlet injection (v36) is still commented out, as in v30.
