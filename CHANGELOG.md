# Changelog

## Unreleased performance branch

Exact stability-only probes for discarded initial CFL attempts, local immutable
parallel-face geometry, and removal of duplicate post-chemistry halo/response
work. Qualified against v0.1.0 with unchanged physical/step policies; a separate
80MPIx3 layout pair also passes. See docs/PERFORMANCE_UPDATE.md. The v0.1.0
tag and bundled baseline methodology artifacts are preserved.

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
