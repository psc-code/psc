
#include <mrc_vec_private.h>

#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <petscvec.h>

#define CE CHKERRABORT(PETSC_COMM_WORLD, ierr)

// ======================================================================
// mrc_vec_petsc subclass

struct mrc_vec_petsc {
  Vec petsc_vec;
  int block_size;
  bool with_petsc_vec;
  bool gotten_petsc_vec;
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
  
  if (!vec->with_array && !sub->with_petsc_vec) {
    //vec->obj.comm
    if (sub->block_size) {
      ierr = VecCreate(petsccomm, &sub->petsc_vec); CE;
      ierr = VecSetSizes(sub->petsc_vec, vec->len, PETSC_DECIDE); CE;
      ierr = VecSetBlockSize(sub->petsc_vec, sub->block_size); CE;
      ierr = VecSetUp(sub->petsc_vec); CE;
    } else {
      ierr = VecCreateMPI(petsccomm, vec->len, PETSC_DETERMINE, &sub->petsc_vec); CE;
    }
  } else if (!sub->petsc_vec) {
    // I'd rather have this in vec_set_sub_set_array, but the calling order doesn't work
    // Petsc needs the len to initialize, and that isn't assigned until mrc_fld_setup.
    // Since 'set array' is called before 'setup' we have to create the petsc vec here instead of there
    assert(vec->arr);
    ierr = VecCreateMPIWithArray(petsccomm, 1, vec->len, PETSC_DECIDE, vec->arr, &sub->petsc_vec); CE;
    // FIXME: I guess our blocksize is 1 for now?
  } else {
    ierr = VecGetBlockSize(sub->petsc_vec, &sub->block_size); CE;
    ierr = VecGetLocalSize(sub->petsc_vec, &vec->len); CE;
  }
}

// ----------------------------------------------------------------------
// mrc_vec_destroy

static void
_mrc_vec_petsc_destroy(struct mrc_vec *vec)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);
  assert(!sub->gotten_petsc_vec);
  int ierr;
  if (!sub->with_petsc_vec) {
  // Petsc won't dealloc memory for WithArray vecs, so it's always safe to call this
    ierr = VecDestroy(&sub->petsc_vec); CE;
  }
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


#define VAR(x) (void *)offsetof(struct mrc_vec_petsc, x)
static struct param mrc_vec_petsc_descr[] = {
  { "block_size"     , VAR(block_size)     , PARAM_INT(0)    },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_vec_petsc_set_petsc_vec
static void
mrc_vec_petsc_set_petsc(struct mrc_vec *vec, Vec invec)
{
  assert(!mrc_vec_petsc(vec)->petsc_vec);
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);
  sub->with_petsc_vec = true;
  sub->petsc_vec = invec;
}

// ----------------------------------------------------------------------
// mrc_vec_petsc_get_petsc_vec
static void
mrc_vec_petsc_get_petsc(struct mrc_vec *vec, Vec *retvec)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);
  assert(!sub->gotten_petsc_vec);
  sub->gotten_petsc_vec = true;
  *retvec = sub->petsc_vec;
}

// ----------------------------------------------------------------------
// mrc_vec_petsc_get_petsc_vec
static void
mrc_vec_petsc_put_petsc(struct mrc_vec *vec, Vec *invec)
{
  struct mrc_vec_petsc *sub = mrc_vec_petsc(vec);
  assert( sub->gotten_petsc_vec && 
	  (sub->petsc_vec == *invec) );
  sub->gotten_petsc_vec = false;
  *invec = NULL;
}

static void
mrc_vec_petsc_axpy(struct mrc_vec *y, double alpha, struct mrc_vec *x)
{
  assert(strcmp(mrc_vec_type(y), mrc_vec_type(x)) == 0);
  int ierr = VecAXPY(mrc_vec_petsc(y)->petsc_vec, 
		     (PetscScalar) alpha, 
		     mrc_vec_petsc(x)->petsc_vec); CE;
}

static void
mrc_vec_petsc_waxpy(struct mrc_vec *w, double alpha, struct mrc_vec *x,  struct mrc_vec *y)
{
  assert(strcmp(mrc_vec_type(y), mrc_vec_type(x)) == 0);
  assert(strcmp(mrc_vec_type(w), mrc_vec_type(x)) == 0);
  int ierr = VecWAXPY(mrc_vec_petsc(w)->petsc_vec, 
		     (PetscScalar) alpha, 
		      mrc_vec_petsc(x)->petsc_vec,
		      mrc_vec_petsc(y)->petsc_vec); CE;
}

static double
mrc_vec_petsc_norm(struct mrc_vec *x)
{
  PetscReal norm;
  int ierr = VecNorm(mrc_vec_petsc(x)->petsc_vec, NORM_2, &norm); CE;
  return (double) norm;
}

static void
mrc_vec_petsc_set(struct mrc_vec *x, double val)
{
  int ierr = VecSet(mrc_vec_petsc(x)->petsc_vec, (PetscScalar) val); CE;
}

static void
mrc_vec_petsc_copy(struct mrc_vec *vec_to, struct mrc_vec *vec_from)
{
  assert(strcmp(mrc_vec_type(vec_to), mrc_vec_type(vec_from)) == 0);
  int ierr = VecCopy(mrc_vec_petsc(vec_from)->petsc_vec,
		     mrc_vec_petsc(vec_to)->petsc_vec); CE;
}


// ----------------------------------------------------------------------
// mrc_vec_petsc_methods
// I think I'm using this right???
// These should probably be renamed as they are very confusing.
// On the other hand, they are palindromic, which I think is super cool.
static struct mrc_obj_method mrc_vec_petsc_methods[] = {
  MRC_OBJ_METHOD("set_petsc_vec", mrc_vec_petsc_set_petsc),
  MRC_OBJ_METHOD("get_petsc_vec", mrc_vec_petsc_get_petsc),
  MRC_OBJ_METHOD("put_petsc_vec", mrc_vec_petsc_put_petsc),
  MRC_OBJ_METHOD("axpy", mrc_vec_petsc_axpy),
  MRC_OBJ_METHOD("waxpy", mrc_vec_petsc_waxpy),
  MRC_OBJ_METHOD("norm", mrc_vec_petsc_norm),
  MRC_OBJ_METHOD("set", mrc_vec_petsc_set),
  MRC_OBJ_METHOD("copy", mrc_vec_petsc_copy),
  {},
};


// ----------------------------------------------------------------------
// mrc_vec subclass "petsc"

struct mrc_vec_ops mrc_vec_petsc_ops = {
  .name                  = "petsc",		
  .methods               = mrc_vec_petsc_methods,
  .size                  = sizeof(struct mrc_vec_petsc),
  .param_descr           = mrc_vec_petsc_descr,
  .create                = mrc_vec_petsc_create,	
  .setup                 = _mrc_vec_petsc_setup,
  .destroy               = _mrc_vec_petsc_destroy,
  .set_array             = mrc_vec_sub_set_array,	
  .get_array             = mrc_vec_sub_get_array,	
  .put_array             = mrc_vec_sub_put_array,
  .axpy                  = mrc_vec_petsc_axpy,
  .waxpy                 = mrc_vec_petsc_waxpy,
  .set                   = mrc_vec_petsc_set,
  .copy                  = mrc_vec_petsc_copy,
};

