# base-case → room5x5: what changed and why it works

Oct 30, 2026 · Krishnaprasad K. A. · Online copy: https://claude.ai/artifact/3yB1zEFuURqthaRovZg6Q1

## Summary

room5x5 builds and runs because it fixes three separate faults that stopped the earlier cases: a Slurm rank layout the current Nek5000 cannot read meshes with, a room grid that did not line up with the outlet grid, and a mixed Intel/gfortran build. It also fixes several silent bugs in the 2022 user code.

The earlier cases are `base-case/` (10 × 10 m room, fluid only) and `base-case2/` (first 5 × 5 m attempt with ppiclF). Neither ran. room5x5 is the 5 × 5 × 3.2 m room rebuilt from the 2022 reference `room_one_v30` on one clean toolchain.

| case | room | what stopped it | logged error |
| --- | --- | --- | --- |
| `base-case/` | 10 × 10 m, 60×19×60 + outlet (68,448 elements) | `--distribution=cyclic:cyclic` over 4–16 nodes | `Error reading .re2 mesh ierr=640` (3 runs) |
| `base-case2/` | 5 × 5 m, 30×19×30 + outlet (17,148 elements) | Room x/z nodes miss the outlet nodes; Intel 2018 Nek linked to gfortran-built ppiclF | `Error invalid tuple sum!` (3 runs); `undefined reference to _gfortran_*` at link |
| `room5x5/` | 5 × 5 m, same element count | All three fixed; see below | Builds clean; mesh check passes |

## Toolchain and build

room5x5 uses one compiler family, pinned source versions, and its own ppiclF build, which removes the link failure seen in base-case2.

| item | base-case / base-case2 | room5x5 |
| --- | --- | --- |
| Compiler + MPI | base-case: `module load gcc openmpi` (unpinned); base-case2: `intel/2018` + `mpif77` + MKL | `gcc/14.2.0` + `openmpi/5.0.10`, `mpif90` |
| Nek5000 | Shared `~/Nek5000` from 2022 | `~/codes/Nek5000`, master commit `d083ef3`, cloned by `setup_nek.sh` |
| ppiclF | Shared `~/ppiclF/source`, built with gfortran | `~/codes/ppiclF_<case>`, commit `72483db`, built per case by `build_ppiclf.sh` |
| Compile flags | `-mcmodel=large -g -traceback`, MKL link group | `-mcmodel=medium`; ppiclF with `-cpp -std=legacy -O2` |
| genbox / genmap | Binaries copied into each case folder | Built once in `~/codes/Nek5000/bin`; `make_mesh.sh` calls them |

Why it matters:

- **Mixed compilers broke the link.** base-case2 compiled Nek5000 with Intel 2018 but linked a ppiclF library compiled with gfortran, so every gfortran runtime call was missing (`_gfortran_st_write`, `_gfortran_transfer_real_write`, and others).
- **ppiclF must be built per case.** ppiclF compiles `ppiclf_user.f` and `PPICLF_USER.h` into `libppiclF.a`. A shared library silently carries another case's particle equations and array sizes.
- **`-std=legacy`** lets gfortran 10+ accept the old argument-type mismatches in ppiclF and Nek user code instead of stopping with errors.
- **Location:** the Nek5000 and ppiclF clones now live under `~/codes/`, away from the 2022 `~/Nek5000` and `~/ppiclF`.

## Mesh and geometry

The room grid in room5x5 now passes through every outlet-box node, so the two boxes share faces and genmap finds 19,295 unique vertices. In base-case2 the boxes did not match, which caused `Error invalid tuple sum!`.

| x/z node (m) | outlet box needs | base-case2 room grid | room5x5 room grid |
| --- | --- | --- | --- |
| 1 | 2.20 | 2.182 (miss) | 2.20 |
| 2 | 2.35 | 2.35 | 2.35 |
| 3 | 2.50 | 2.518 (miss) | 2.50 |
| 4 | 2.65 | 2.686 (miss) | 2.65 |
| 5 | 2.80 | 2.854 (miss) | 2.80 |

- **base-case2 grid:** it kept the 10 × 10 room spacing of 0.1679 m and simply cut it off at 5 m, ending with uneven 4.70 / 4.85 / 5.00 m nodes.
- **room5x5 grid:** uniform 0.1692 m away from the outlet, 0.15 m across the 2.2–2.8 m outlet band in x and z, and 30×19×30 + 4×3×4 = 17,148 elements.
- **SIZE:** `lelg = 17148`, `lpmin = 32`, `lx1 = 6`. base-case used `lelg = 68448`, `lpmin = 64` for the 10 × 10 room.
- **Check:** `make_mesh.sh` builds `room.re2` and `room.ma2` and prints the vertex count. Any number above 19,295 means the boxes no longer line up.
- **Geometry in one place:** room size, outlet and inlet slot positions sit in one block at the top of `usrdat2`. Ceiling elements under the outlet are found by coordinates, not hard-coded element numbers, so the mesh can change without editing the boundary-condition code.

## Physics and code fixes

The physics model is unchanged from the 2022 code: dynamic Smagorinsky LES, an outlet sponge, Langevin sub-grid particle velocity, and Stokes settling. The fixes below correct bugs that the base cases inherited from that code. Several were silent: Intel tolerated them, but gfortran either crashed or gave wrong results.

| area | 2022 code (base cases, v30) | room5x5 | effect |
| --- | --- | --- | --- |
| Inlet speed | Fixed at 7.1428 m/s, correct only for the 10 × 10 room | `uin = ACH·V/3600 / (4·A_inlet·sin θ)` from `userParam06` (ACH) and `userParam07` (vane angle) | Gives 1.79 m/s for the 5 × 5 room at ACH 5 and 40°, the same airflow per room volume |
| `userbc` | `ux` kept its value from the previous call | `ux` zeroed every call | Removes cross-talk between inlet faces |
| Random numbers | `ran2` / `ppiclf_ran` state not SAVEd; Langevin generator never seeded | State SAVEd; Langevin generator seeded per rank | Undefined values under gfortran; before, every rank drew the same particle noise |
| Element loops | Loops ran to `lelv` (array size) | Loops run to `nelv` (elements actually on the rank) | No reads of unset elements |
| Wall functions | Outlet-box elements indexed past array bounds | Outlet-box elements skipped | No out-of-bounds writes |
| `userchk` | `gop` used Nek's mesh-velocity arrays as scratch; wall-distance arrays not SAVEd; a declaration missed a comma | Own work array, SAVEd arrays, declaration fixed | Mesh-velocity arrays and wall distances stay intact |
| Particles switch | ppiclF always on | `userParam05` = 0 skips ppiclF entirely | Fluid-only spin-up runs without particle memory or cost |

Inlet injection (the v36 style) is still commented out, as in v30.

## Run configuration

The Slurm script no longer uses cyclic rank placement, which is the single change that clears base-case's `Error reading .re2 mesh ierr=640`.

| item | base-case / base-case2 | room5x5 |
| --- | --- | --- |
| Rank placement | `--distribution=cyclic:cyclic` | Default (block) placement |
| Layout | base-case: 64 ranks over 4 nodes in the logged runs; base-case2: 64 ranks | 1 node, 64 ranks (32 also allowed) |
| QOS / time | `bala1s` 10 min, or `bala1s-b` | `bala1s-b`, 4 days; hpg-dev test uses `bala1s`, 30 min |
| MPI launch | Custom UCX settings (`OMPI_MCA_pml=ucx`, `UCX_TLS`) | `srun --mpi=${HPC_PMIX:-pmix}` with module defaults |
| Stop condition | `numSteps = 10000000` | `stopAt = endTime` (600 s; 100 s in room5x5_v2) |
| Run controls | One `userparam01` | `userParam01`–`10`: particle output, particles on/off, ACH, vane angle, particle count, re-seeding |

Why cyclic placement fails: the current Nek5000 mesh reader uses one aggregator rank per node and assumes each node holds a contiguous block of ranks. Cyclic placement interleaves ranks across nodes. The logs show this pattern, with tasks 0, 4, 8, … on one node and 1, 5, 9, … on the next. The read then returns error 10 on all 64 ranks: 64 × 10 = 640.

Particle memory is now checked at start-up. `PPICLF_LPART = 160000` per rank, and `room.usr` stops if particles per rank exceed LPART/2. So 5 M particles need 64 ranks.

## Status and open items

room5x5 is proven to build, not yet proven to run at production size: no lx1 = 6 run has finished yet.

Tested so far:

- **Cloud workspace, reduced size** (lx1 = 4, 4 ranks): fluid-only spin-up and a 40k-particle restart both ran with no NaNs.
- **HiPerGator, production size, 30 Oct 2026:** mesh check passed (19,295 vertices, 17,148 elements); ppiclF and `nek5000` compiled with gcc 14.2.

Jobs queued on 30 Oct 2026:

| job | case | purpose | settings | state |
| --- | --- | --- | --- | --- |
| 48367635 | `room5x5/` | hpg-dev smoke test | 200 steps, 32 ranks, 40° vanes | Pending, estimated start 22:41 |
| 48367310 | `room5x5_v2/` | Fluid spin-up to t = 100 s | 64 ranks, bala1s-b, 60° vanes (inlet about 1.55 m/s) | Pending |

What a healthy log shows: the `room5x5: ACH … vane angle … -> inlet speed …` line, steps advancing, and CFL `C=` below about 0.6.

Open items:

- [ ] After job 48367635 ends, restore `room5x5/room.par` from `room.par.prod` (it is set to 200 steps for the test).
- [ ] Commit `SIZE` and `room.par` for the room5x5 cases; `.gitignore` currently excludes them, so a fresh clone cannot build.
- [ ] Point `room5x5_particles/` scripts at `~/codes/` like the other two cases.
- [ ] Decide the ACH for the 5 × 5 room: at ACH 5 the inlet Reynolds number is about 4× lower than in the 10 × 10 room.
- [ ] Decide on turning inlet injection back on (v36 style) and on periodic boundaries for the multi-room stages.
