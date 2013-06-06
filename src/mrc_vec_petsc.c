
#include <mrc_vec_private.h>

#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <petscvec.h>


// ----------------------------------------------------------------------
// mrc_vec_setup

static void
_mrc_vec_setup(struct mrc_vec *vec)
{
  if (!vec->with_array) {
    vec->arr = calloc(vec->len, vec->size_of_type);
  }
}

// ----------------------------------------------------------------------
// mrc_vec_destroy

static void
_mrc_vec_destroy(struct mrc_vec *vec)
{
  if (!vec->with_array) {
    free(vec->arr);
  }
  vec->arr = NULL;
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
  return vec->arr;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_put_array

static void
mrc_vec_sub_put_array(struct mrc_vec *vec, void *arr)
{
  assert(arr == vec->arr);
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
  .set_array             = mrc_vec_sub_set_array,	
  .get_array             = mrc_vec_sub_get_array,	
  .put_array             = mrc_vec_sub_put_array,	
};

