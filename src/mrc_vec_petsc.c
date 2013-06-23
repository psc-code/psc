
#include <mrc_vec_private.h>

#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <petscvec.h>

#define CE CHKERRABORT(vec->obj.comm, ierr)

// ======================================================================
// mrc_vec_petsc subclass

struct mrc_vec_petsc {
  Vec petsc_vec;
};

#define mrc_vec_petsc(vec) mrc_to_subobj(vec, struct mrc_vec_petsc)

// ----------------------------------------------------------------------
// mrc_vec_setup

static void
_mrc_vec_petsc_setup(struct mrc_vec *vec)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);

  int ierr;
  MPI_Comm petsccomm;
  // This seems to be a valid way to handle the MPI_COMM_NULL edge case,
  // but I'm not entirely sure how VecCreateMPI behaves on PETSC_COMM_SELF
  if (vec->obj.comm != MPI_COMM_NULL) {
    petsccomm = vec->obj.comm;
  } else {
    petsccomm = PETSC_COMM_SELF;
  }
  
  if (!vec->with_array) {
    //vec->obj.comm
    ierr = VecCreateMPI(petsccomm, vec->len, PETSC_DETERMINE, &sub->petsc_vec); CE;
  } else {
    // I'd rather have this in vec_set_sub_set_array, but the calling order doesn't work
    // Petsc needs the len to initialize, and that isn't assigned until mrc_fld_setup.
    // Since 'set array' is called before 'setup' we have to create the petsc vec here instead of there
    assert(vec->arr);
    ierr = VecCreateMPIWithArray(petsccomm, 1, vec->len, PETSC_DECIDE, vec->arr, &sub->petsc_vec); CE;
    // FIXME: I guess our blocksize is 1 for now?
  }
}

// ----------------------------------------------------------------------
// mrc_vec_destroy

static void
_mrc_vec_petsc_destroy(struct mrc_vec *vec)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);

  int ierr;
  // Petsc won't dealloc memory for WithArray vecs, so it's always safe to call this
  ierr = VecDestroy(&sub->petsc_vec); CE;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_set_array

static void
mrc_vec_sub_set_array(struct mrc_vec *vec, void *arr)
{
  assert(!vec->arr);  
  vec->arr = arr;
  vec->with_array = true;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_get_array

static void *
mrc_vec_sub_get_array(struct mrc_vec *vec)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);
  int ierr;
  ierr = VecGetArray(sub->petsc_vec, (PetscScalar **) &vec->arr); CE;
  return vec->arr;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_put_array

static void
mrc_vec_sub_put_array(struct mrc_vec *vec, void *arr)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);
  int ierr;
  assert(arr == vec->arr);
  ierr = VecRestoreArray(sub->petsc_vec, (PetscScalar **) &vec->arr); CE;
}

// ----------------------------------------------------------------------
// mrc_vec_petsc_create

static void
mrc_vec_petsc_create(struct mrc_vec *vec)
{
  // Set type from petsc
  vec->size_of_type = sizeof(PetscScalar);
}

// ----------------------------------------------------------------------
// mrc_vec subclass "petsc"

struct mrc_vec_ops mrc_vec_petsc_ops = {
  .name                  = "petsc",		
  .size                  = sizeof(struct mrc_vec_petsc),
  .create                = mrc_vec_petsc_create,	
  .setup                 = _mrc_vec_petsc_setup,
  .destroy               = _mrc_vec_petsc_destroy,
  .set_array             = mrc_vec_sub_set_array,	
  .get_array             = mrc_vec_sub_get_array,	
  .put_array             = mrc_vec_sub_put_array,	
};

