
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

// ----------------------------------------------------------------------
// mrc_vec_setup

static void
_mrc_vec_petsc_setup(struct mrc_vec *vec)
{
  int ierr;
  // FIXME: if with_array is given we should do a VecCreateMPIWithArray
  if (!vec->with_array) {
    //     vec->arr = calloc(vec->len, vec->size_of_type);
    ierr = VecCreateMPI(vec->obj.comm, vec->len, PETSC_DETERMINE, (Vec *) &vec->_priv_arr); CE;
  }
}

// ----------------------------------------------------------------------
// mrc_vec_destroy

static void
_mrc_vec_petsc_destroy(struct mrc_vec *vec)
{
  int ierr;
  if (!vec->with_array) {
    ierr = VecDestroy((Vec *) &vec->_priv_arr); CE;
  }
  vec->_priv_arr = NULL;
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
  int ierr;
  ierr = VecGetArray((Vec) vec->_priv_arr, (PetscScalar **) &vec->arr); CE;
  return vec->arr;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_put_array

static void
mrc_vec_sub_put_array(struct mrc_vec *vec, void *arr)
{
  int ierr;
  assert(arr == vec->arr);
  ierr = VecRestoreArray((Vec) vec->_priv_arr, (PetscScalar **) &vec->arr); CE;
}

// ======================================================================
// mrc_vec_petsc subclass



// ----------------------------------------------------------------------
// mrc_vec_petsc_create

static void
mrc_vec_petsc_create(struct mrc_vec *vec)
{
  // Set type from petsc
  vec->size_of_type = sizeof(PetscScalar);
}

struct mrc_vec_ops mrc_vec_petsc_ops = {
  .name                  = "petsc",			
  .create                = mrc_vec_petsc_create,	
  .setup                 = _mrc_vec_petsc_setup,
  .destroy               = _mrc_vec_petsc_destroy,
  .set_array             = mrc_vec_sub_set_array,	
  .get_array             = mrc_vec_sub_get_array,	
  .put_array             = mrc_vec_sub_put_array,	
};

