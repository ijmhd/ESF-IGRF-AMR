# Changelog

## main audit fixes - 2026-11-16

* Bind every newly written native AMR checkpoint to the complete numerical
  source ID with an `ESFSRC01` envelope covered by its checksum, including
  non-night and legacy-mask AMR paths. Unbound older native/static imports must
  use their original executable; no automatic conversion or bypass is supplied.
* Use the shared saved-frame reader for 2D slices and movies, including missing
  two-level masks, incomplete CSV tails and O+/electron-specific file needs;
  refuse an existing figure/provenance prefix before writing.
* Separate segment failure retries from restored cumulative history while
  retaining the original checkpoint/completion totals.
* Record the targeted short-case gates in `validation/audit_fixes/`, update the
  Chinese tutorial and publish an explicit old-checkpoint compatibility table.
  Physics, evolution algorithms, CLI defaults, example inputs and build defaults
  remain unchanged. The `v0.1.0` tag, release attachments and original validation
  records remain frozen; no new tag or replacement historical result is created.

## v0.1.0 - 2026-11-16

First frozen research/teaching release of the independent C11 ESF solver.

* Dipole/IGRF-14 traced conjugate geometry, O+ continuity/parallel momentum,
  integrated full Pedersen/Hall potential and constrained PETSc AMR backend.
* Two/three physical transverse levels, segmented density refinement,
  solution-based wall/tip/branch criteria, conservative transfer, proper nesting,
  physical-time coarsening hysteresis and optional early-mode protection.
* Optional reduced transported-night NO+/O2+ chemistry and common electron
  pressure, with prescribed backgrounds and explicit physical limitations.
* Owned molecular migration and stage-valid caches, cache-friendly independent
  OpenMP kernels and separate attempted-step/CFL counters. The exact source and
  completed performance evidence are recorded in `validation/`.
* Adjacent responsibility comments for every runtime C function, data/ownership
  contracts, all79 CLI controls, explicit small cases, MPI/Slurm/restart/plot
  instructions, Chinese learning guide, LaTeX methods and full source atlas.

The initial annotation/formatting pass was token-identical to its pre-annotation
source. Subsequent intentional traversal/status/cache changes have their own
source identities and developed-state comparisons. The released source is
identified by `validation/source_identity.json`; a Git commit additionally binds
the documentation and test sources. The source-file manifest excludes itself
to avoid a circular hash.

This release does not contain general octree AMR, molecular AMR, a new MPI tile
decomposition or a change of the retained CFL/solver/refinement policies. Its
timing qualification is a150s continuation of a developed state, not an event
calibration or a full strong-scaling study. See the validation and audit guides.
