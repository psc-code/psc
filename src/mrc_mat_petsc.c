
#include "mrc_mat_private.h"

#include <petscmat.h>

#include <stdlib.h>

#define CE CHKERRABORT(PETSC_COMM_WORLD, ierr)

// ======================================================================
// mrc_mat "petsc"

struct mrc_mat_petsc {
  Mat mat;
};

#define mrc_mat_petsc(mat) mrc_to_subobj(mat, struct mrc_mat_petsc)

// ----------------------------------------------------------------------
// mrc_mat_petsc_create

static void
mrc_mat_petsc_create(struct mrc_mat *mat)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  ierr = MatCreate(mrc_mat_comm(mat), &sub->mat); CE;
}

// ----------------------------------------------------------------------
// mrc_mat_petsc_setup

static void
mrc_mat_petsc_setup(struct mrc_mat *mat)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  mprintf("set sizes %d %d\n", mat->m, mat->n);
  ierr = MatSetSizes(sub->mat, mat->m, mat->n, PETSC_DECIDE, PETSC_DECIDE); CE;
  ierr = MatSetUp(sub->mat); CE;
}

// ----------------------------------------------------------------------
// mrc_mat_petsc_destroy

static void
mrc_mat_petsc_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  ierr = MatDestroy(&sub->mat); CE;
}

// ----------------------------------------------------------------------
// mrc_mat_petsc_add_value

static void
mrc_mat_petsc_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  ierr = MatSetValue(sub->mat, row_idx, col_idx, val, ADD_VALUES); CE;
}

// ----------------------------------------------------------------------
// mrc_mat_petsc_assemble

static void
mrc_mat_petsc_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  ierr = MatAssemblyBegin(sub->mat, MAT_FINAL_ASSEMBLY); CE;
  ierr = MatAssemblyEnd(sub->mat, MAT_FINAL_ASSEMBLY); CE;
}

// ----------------------------------------------------------------------
// mrc_mat_petsc_apply

typedef Vec (*fgp_t)(struct mrc_fld *);
typedef void (*fpp_t)(struct mrc_fld *, Vec *);

static void
mrc_mat_petsc_apply(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  fgp_t fld_get_petsc = (fgp_t) mrc_fld_get_method(x, "get_petsc_vec");
  fpp_t fld_put_petsc = (fpp_t) mrc_fld_get_method(x, "put_petsc_vec");

  Vec xvec = fld_get_petsc(x);
  Vec yvec = fld_get_petsc(y);

  ierr = MatMult(sub->mat, xvec, yvec); CE;

  fld_put_petsc(x, &xvec);
  fld_put_petsc(y, &yvec);
}

// ----------------------------------------------------------------------
// mrc_mat_petsc_apply_in_place

static void
mrc_mat_petsc_apply_in_place(struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_petsc *sub = mrc_mat_petsc(mat);
  int ierr;

  fgp_t fld_get_petsc = (fgp_t) mrc_fld_get_method(x, "get_petsc_vec");
  fpp_t fld_put_petsc = (fpp_t) mrc_fld_get_method(x, "put_petsc_vec");

  Vec xvec, yvec;
  xvec = fld_get_petsc(x);
  ierr = VecDuplicate(xvec, &yvec); CE;

  ierr = MatMult(sub->mat, xvec, yvec); CE;
  ierr = VecCopy(yvec, xvec); CE;

  ierr = VecDestroy(&yvec); CE;
  fld_put_petsc(x, &xvec);
}

// ----------------------------------------------------------------------
// mrc_mat subclass "petsc"

struct mrc_mat_ops mrc_mat_petsc_ops = {
  .name                  = "petsc",		
  .size                  = sizeof(struct mrc_mat_petsc),
  .create                = mrc_mat_petsc_create,
  .setup                 = mrc_mat_petsc_setup,
  .destroy               = mrc_mat_petsc_destroy,
  .add_value             = mrc_mat_petsc_add_value,
  .assemble              = mrc_mat_petsc_assemble,
  .apply                 = mrc_mat_petsc_apply,
  .apply_in_place        = mrc_mat_petsc_apply_in_place,
};

