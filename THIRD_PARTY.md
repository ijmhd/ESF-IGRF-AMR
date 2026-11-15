# Provenance and dependencies

The independent C implementation is GPL-3.0-only; see LICENSE. No TELF, Athena++
or SWMF source is bundled. Their papers inform the physical/numerical comparison,
not an assertion that this is their solver. The annotated source and methods
state exactly which approximations are implemented.

## IGRF

Official NOAA/IAGA IGRF-14 coefficients:
https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt

SHA256: 8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c
Retrieved2026-11-10. The C harmonic evaluator is independently derived. The
coefficient provenance is distinct from the software licenses.

## PDM limiter

The limiter follows the previously audited GPLv3 PDM4 adaptation and PDMU7.m in
https://github.com/bzhang1982/gamera_matlab . The mapped physical-volume moments
and four-cell cubic weights are implemented for this solver. Current PDM4 is
not PDM7 and does not imply a fourth-order coupled time method.

## Reaction and collision formulas

The thermal low-effective-temperature reaction fits are traced to GEMINI3D
commit d70446a712d47e7617b25c033bfcace64207590f,
https://github.com/gemini3d/gemini3d/blob/d70446a712d47e7617b25c033bfcace64207590f/src/sources/sources.f90
(source blob a96d9ac3e6913a491eeb72a7e5a3140158d62fa0, Apache-2.0 project).
The C local positivity-preserving reaction composition is reduced and independent;
this bundle contains no GEMINI Fortran source and is not a GEMINI model.

Molecular collision/conductivity expression provenance:
https://github.com/gemini3d/glow/blob/master/conduct.f90
The exact retained expressions and empirical electron multiplier4 are printed
in the frozen methodology and source. No GLOW program is bundled. Prescribed
O+ collision input policies are separate; a synthetic teaching CSV is not MSIS.

## External software

PETSc3.26.1 real-double is the qualified optional backend, under its BSD license:
https://petsc.org/release/install/license/ . PETSc/MPI/BLAS/LAPACK/OpenMP runtimes
are not bundled and retain their own licenses. Qualification build flags and
version/resource evidence are in validation/qualification.json. Python standard
library runs the launcher; optional NumPy/matplotlib readers are external.
LaTeX packages only build the methods documents.

The bundle includes no observed ISR data, TELF HDF5 data, cluster credentials,
large historical checkpoints or native geometry caches. Related-work references
are in methodology/references.bib. Measured software results are this project's
independent qualification, with its explicit scope and limitations.
