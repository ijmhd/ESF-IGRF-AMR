/*
 * Distributed PETSc backend for the uniform-grid reference solver; the AMR solver has its own
 * constrained-mesh assembly.
 * See docs/METHODS.md and the methodology/ LaTeX project for equations and call order.
 */
#include "solver.h"
#include <petscksp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    Mat a;
    Vec b, x, residual, all;
    VecScatter scatter;
    KSP ksp;
    PetscInt first, last;
    int solves, last_iterations;
} ParallelPotential;

/* PETSc owns distributed rows, Krylov vectors and multigrid operations.
 * The existing Q1 assembly supplies exactly the same owned matrix rows,
 * including Hall, periodic seam and regional boundary terms. */
/**
 * Set a PETSc backend default only when that option has not already been explicitly set.
 */
static PetscErrorCode default_option(const char *name, const char *value) {
    PetscBool present;
    PetscCall(PetscOptionsHasName(NULL, NULL, name, &present));
    if (!present)
        PetscCall(PetscOptionsSetValue(NULL, name, value));
    return 0;
}
/**
 * Allocate distributed uniform-grid PETSc matrices/vectors/scatter and configure the reference
 * potential solve.
 */
int solver3_petsc_create(Solver3 *s) {
    ParallelPotential *c = alloc3(s, 1, sizeof(*c));
    s->petsc_context = c;
    PC pc;
    PetscCall(MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, s->p.n, s->p.n, 9, NULL, 9,
                           NULL, &c->a));
    PetscCall(MatGetOwnershipRange(c->a, &c->first, &c->last));
    s->p.first_row = c->first;
    s->p.last_row = c->last;
    PetscCall(MatSetOption(c->a, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
    PetscCall(MatCreateVecs(c->a, &c->x, &c->b));
    PetscCall(VecDuplicate(c->b, &c->residual));
    PetscCall(VecSet(c->x, 0));
    PetscCall(VecScatterCreateToAll(c->x, &c->scatter, &c->all));
    PetscCall(KSPCreate(MPI_COMM_WORLD, &c->ksp));
    PetscCall(KSPSetOptionsPrefix(c->ksp, "esf_"));
    PetscCall(KSPSetOperators(c->ksp, c->a, c->a));
    PetscCall(KSPSetType(c->ksp, KSPFGMRES));
    PetscCall(KSPGMRESSetRestart(c->ksp, 40));
    PetscCall(KSPSetPCSide(c->ksp, PC_RIGHT));
    PetscCall(KSPSetInitialGuessNonzero(c->ksp, PETSC_TRUE));
    PetscCall(KSPGetPC(c->ksp, &pc));
    PetscCall(PCSetType(pc, PCGAMG));
    PetscCall(PCGAMGSetNSmooths(pc, 1));
    PetscCall(default_option("-esf_mg_levels_ksp_type", "richardson"));
    PetscCall(default_option("-esf_mg_levels_ksp_max_it", "2"));
    PetscCall(default_option("-esf_mg_levels_pc_type", "sor"));
    PetscCall(default_option("-esf_pc_gamg_esteig_ksp_type", "gmres"));
    PetscCall(PCGAMGSetReuseInterpolation(pc, PETSC_TRUE));
    PetscCall(KSPSetFromOptions(c->ksp));
    /* Always test the unpreconditioned system. In particular, a small
     * preconditioned norm is not accepted in place of the physical residual. */
    PetscCall(KSPSetNormType(c->ksp, KSP_NORM_UNPRECONDITIONED));
    return 0;
}

/**
 * Assemble current uniform-grid coefficients, apply boundary conditions, solve and verify the
 * true residual before distributing potential.
 */
int solver3_petsc_solve(Solver3 *s) {
    ParallelPotential *c = s->petsc_context;
    Potential *p = &s->p;
    double tick = wall3();
    int columns = s->c.grid.nx * s->c.grid.ny;
    for (int i = 0; i < columns; i++) {
        s->coeff[6 * i + 4] -= s->coeff[6 * i] * s->background_potential_gradient;
        s->coeff[6 * i + 5] -= s->coeff[6 * i + 2] * s->background_potential_gradient;
    }
    potential_assemble(p, s->coeff);
    PetscScalar *b;
    PetscCall(VecGetArray(c->b, &b));
    for (PetscInt row = c->first; row < c->last; row++) {
        PetscInt indices[9], nnz = p->row[row + 1] - p->row[row];
        for (PetscInt k = 0; k < nnz; k++)
            indices[k] = p->col[p->row[row] + k];
        PetscCall(MatSetValues(c->a, 1, &row, nnz, indices, p->a + p->row[row], INSERT_VALUES));
        b[row - c->first] = p->rhs[row];
    }
    PetscCall(VecRestoreArray(c->b, &b));
    PetscCall(MatAssemblyBegin(c->a, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(c->a, MAT_FINAL_ASSEMBLY));
    s->timing[PROF_ASSEMBLY] += wall3() - tick;
    tick = wall3();
    PetscReal scale;
    PetscCall(VecNorm(c->b, NORM_INFINITY, &scale));
    p->iterations = 0;
    p->residual = 0;
    if (scale == 0) {
        PetscCall(VecSet(c->x, 0));
    } else {
        /* Slowly changing coefficients permit reuse of a preconditioner.
         * Rebuild regularly, after expensive solves, or on a failed true
         * residual check. The matrix and residual tolerance never freeze. */
        PetscBool reuse = (c->solves % 32 && c->last_iterations < 20) ? PETSC_TRUE : PETSC_FALSE;
        for (int attempt = 0; attempt < 2; attempt++) {
            PetscCall(KSPSetReusePreconditioner(c->ksp, attempt ? PETSC_FALSE : reuse));
            PetscCall(KSPSetTolerances(c->ksp, 1e-14, .25 * s->c.rtol * scale, PETSC_DEFAULT, 400));
            PetscCall(KSPSolve(c->ksp, c->b, c->x));
            PetscInt iterations;
            PetscCall(KSPGetIterationNumber(c->ksp, &iterations));
            p->iterations += iterations;
            PetscCall(MatMult(c->a, c->x, c->residual));
            PetscCall(VecAYPX(c->residual, -1, c->b));
            PetscReal norm;
            PetscCall(VecNorm(c->residual, NORM_INFINITY, &norm));
            p->residual = norm / scale;
            if (isfinite(p->residual) && p->residual <= 1.05 * s->c.rtol)
                break;
        }
        if (!isfinite(p->residual) || p->residual > 1.05 * s->c.rtol)
            return -1;
    }
    c->last_iterations = p->iterations;
    c->solves++;
    s->potential_iterations += p->iterations;
    s->max_potential_residual = fmax(s->max_potential_residual, p->residual);
    s->timing[PROF_SOLVE] += wall3() - tick;
    tick = wall3();
    PetscCall(VecScatterBegin(c->scatter, c->x, c->all, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(c->scatter, c->x, c->all, INSERT_VALUES, SCATTER_FORWARD));
    const PetscScalar *x;
    PetscCall(VecGetArrayRead(c->all, &x));
    for (int i = 0; i < p->n; i++)
        s->phi[i] = x[i] + s->background_potential_gradient *
                               ((double)(i % (s->c.grid.nx + 1)) / s->c.grid.nx - .5);
    PetscCall(VecRestoreArrayRead(c->all, &x));
    s->timing[PROF_BCAST] += wall3() - tick;
    return 0;
}

/**
 * Warm-start distributed uniform potential DOFs from the current canonical phi.
 */
int solver3_petsc_seed(Solver3 *s) {
    ParallelPotential *c = s->petsc_context;
    PetscScalar *x;
    PetscCall(VecGetArray(c->x, &x));
    for (PetscInt row = c->first; row < c->last; row++)
        x[row - c->first] =
            s->phi[row] - s->background_potential_gradient *
                              ((double)(row % (s->c.grid.nx + 1)) / s->c.grid.nx - .5);
    PetscCall(VecRestoreArray(c->x, &x));
    return 0;
}

/**
 * Destroy all PETSc objects associated with the uniform-grid backend.
 */
void solver3_petsc_destroy(Solver3 *s) {
    ParallelPotential *c = s->petsc_context;
    PetscCallAbort(MPI_COMM_WORLD, VecScatterDestroy(&c->scatter));
    PetscCallAbort(MPI_COMM_WORLD, VecDestroy(&c->all));
    PetscCallAbort(MPI_COMM_WORLD, VecDestroy(&c->x));
    PetscCallAbort(MPI_COMM_WORLD, VecDestroy(&c->b));
    PetscCallAbort(MPI_COMM_WORLD, VecDestroy(&c->residual));
    PetscCallAbort(MPI_COMM_WORLD, KSPDestroy(&c->ksp));
    PetscCallAbort(MPI_COMM_WORLD, MatDestroy(&c->a));
    free(c);
    s->petsc_context = NULL;
}
