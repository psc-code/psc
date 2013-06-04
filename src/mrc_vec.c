
#include <mrc_vec_private.h>

#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_vec

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
// mrc_vec_set_array
//
// instead of allocating memory, use the pointer we provide

void
mrc_vec_set_array(struct mrc_vec *vec, void *arr)
{
  assert(!vec->arr);
  vec->arr = arr;
  vec->with_array = true;
}

// ----------------------------------------------------------------------
// mrc_vec_get_array
//
// gets the pointer to the allocated storage in the vector

void *
mrc_vec_get_array(struct mrc_vec *vec)
{
  assert(mrc_vec_is_setup(vec));

  return vec->arr;
}

// ----------------------------------------------------------------------
// mrc_vec_put_array
//
// indicates that we're done using the data pointer obtained from
// mrc_vec_get_array

void
mrc_vec_put_array(struct mrc_vec *vec, void *arr)
{
  assert(arr == vec->arr);
}

// ======================================================================
// mrc_vec subclasses

#define MAKE_MRC_VEC_TYPE(type, TYPE)			\
							\
  static void						\
  mrc_vec_##type##_create(struct mrc_vec *vec)		\
  {							\
    vec->size_of_type = sizeof(type);			\
  }							\
  							\
  static struct mrc_vec_ops mrc_vec_##type##_ops = {	\
    .name                  = #type,			\
    .create                = mrc_vec_##type##_create,	\
  };							\

MAKE_MRC_VEC_TYPE(float, FLOAT)
MAKE_MRC_VEC_TYPE(double, DOUBLE)
MAKE_MRC_VEC_TYPE(int, INT)

// ----------------------------------------------------------------------
// mrc_vec_init

static void
mrc_vec_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_vec, &mrc_vec_float_ops);
  mrc_class_register_subclass(&mrc_class_mrc_vec, &mrc_vec_double_ops);
  mrc_class_register_subclass(&mrc_class_mrc_vec, &mrc_vec_int_ops);
}

// ----------------------------------------------------------------------
// mrc_class_mrc_vec

#define VAR(x) (void *)offsetof(struct mrc_vec, x)
static struct param mrc_vec_descr[] = {
  { "len"             , VAR(len)          , PARAM_INT(0)          },

  { "size_of_type"    , VAR(size_of_type) , MRC_VAR_INT           },
  { "with_array"      , VAR(with_array)   , MRC_VAR_BOOL          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_vec class description

struct mrc_class_mrc_vec mrc_class_mrc_vec = {
  .name         = "mrc_vec",
  .size         = sizeof(struct mrc_vec),
  .param_descr  = mrc_vec_descr,
  .init         = mrc_vec_init,
  .setup        = _mrc_vec_setup,
  .destroy      = _mrc_vec_destroy,
};

