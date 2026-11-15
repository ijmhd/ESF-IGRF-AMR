#!/bin/bash
# Example only: set the local partition/account and load your compiler/MPI/PETSc.
#SBATCH --job-name=esf-teaching
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1500M
#SBATCH --time=00:30:00
#SBATCH --output=esf-%j.log
set -euo pipefail
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_DYNAMIC=false OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Configure MPI rank domains with at least OMP_NUM_THREADS CPUs per rank.
# MPI placement flags differ between vendors; see docs/RUNNING.md.
python3 tools/run_case.py hk-amr3 --binary build-omp/esf3d \
  --ranks "${SLURM_NTASKS}" --threads "${OMP_NUM_THREADS}" \
  --out "runs/slurm-${SLURM_JOB_ID}"
