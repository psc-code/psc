
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

#define mrc_vec_ops(vec) ((struct mrc_vec_ops *) vec->obj.ops)

// Public wrappers for math/data ops.
void 
mrc_vec_axpy(struct mrc_vec *y, double alpha, struct mrc_vec *x)
{
  assert(mrc_vec_ops(y)->axpy);
  mrc_vec_ops(y)->axpy(y, alpha, x);
}

void
mrc_vec_waxpy(struct mrc_vec *w, double alpha, struct mrc_vec *x, struct mrc_vec *y)
{
  assert(mrc_vec_ops(w)->waxpy);
  mrc_vec_ops(w)->waxpy(w, alpha, x, y);
}

void 
mrc_vec_axpby(struct mrc_vec *y, double alpha, struct mrc_vec *x, double beta)
{
  assert(mrc_vec_ops(y)->axpby);
  mrc_vec_ops(y)->axpby(y, alpha, x, beta);
}

void 
mrc_vec_set(struct mrc_vec *x, double alpha)
{
  assert(mrc_vec_ops(x)->set);
  mrc_vec_ops(x)->set(x, alpha);
}

void
mrc_vec_copy(struct mrc_vec *vec_to, struct mrc_vec *vec_from)
{
  assert(mrc_vec_ops(vec_to)->copy);
  mrc_vec_ops(vec_to)->copy(vec_to, vec_from);
}


// ----------------------------------------------------------------------
// mrc_vec_setup

static void
_mrc_vec_setup(struct mrc_vec *vec)
{
  if (!vec->with_array) {
    vec->arr = calloc(vec->len, vec->size_of_type);

#ifdef MRC_VEC_INIT_NAN
    memset(vec->arr, 0xff, vec->len * vec->size_of_type);
#endif
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
  struct mrc_vec_ops *ops = mrc_vec_ops(vec);
  assert(ops && ops->set_array);
  ops->set_array(vec, arr);
}

// ----------------------------------------------------------------------
// mrc_vec_replace_array
//
// replace our own allocated memory by the pointer provided (and free our mem)

void
mrc_vec_replace_array(struct mrc_vec *vec, void *arr)
{
  struct mrc_vec_ops *ops = mrc_vec_ops(vec);
  assert(ops && ops->replace_array);
  ops->replace_array(vec, arr);
}

// ----------------------------------------------------------------------
// mrc_vec_get_array
//
// gets the pointer to the allocated storage in the vector

void *
mrc_vec_get_array(struct mrc_vec *vec)
{
  assert(mrc_vec_is_setup(vec));

  struct mrc_vec_ops *ops = mrc_vec_ops(vec);
  assert(ops && ops->get_array);
  return ops->get_array(vec);
}

// ----------------------------------------------------------------------
// mrc_vec_put_array
//
// indicates that we're done using the data pointer obtained from
// mrc_vec_get_array

void
mrc_vec_put_array(struct mrc_vec *vec, void *arr)
{
  assert(mrc_vec_is_setup(vec));

  struct mrc_vec_ops *ops = mrc_vec_ops(vec);
  assert(ops && ops->put_array);
  return ops->put_array(vec, arr);
}

// ----------------------------------------------------------------------
// mrc_vec_sub_set_array

static void
mrc_vec_sub_set_array(struct mrc_vec *vec, void *arr)
{
  assert(arr);
  assert(!vec->arr);
  vec->arr = arr;
  vec->with_array = true;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_replace_array

static void
mrc_vec_sub_replace_array(struct mrc_vec *vec, void *arr)
{
  assert(arr);
  free(vec->arr);
  vec->arr = arr;
}

// ----------------------------------------------------------------------
// mrc_vec_sub__array

static void *
mrc_vec_sub_get_array(struct mrc_vec *vec)
{
  assert(mrc_vec_is_setup(vec));
  return vec->arr;
}

// ----------------------------------------------------------------------
// mrc_vec_sub_put_array

static void
mrc_vec_sub_put_array(struct mrc_vec *vec, void *arr)
{
  assert(mrc_vec_is_setup(vec));
  assert(arr == vec->arr);
}

// ----------------------------------------------------------------------
// mrc_vec_len

int
mrc_vec_len(struct mrc_vec *x)
{
  return x->len;
}

// ----------------------------------------------------------------------
// mrc_vec_size_of_type

int
mrc_vec_size_of_type(struct mrc_vec *x)
{
  return x->size_of_type;
}

// ======================================================================
// mrc_vec subclasses

#define MAKE_MRC_VEC_TYPE(type, TYPE)					\
									\
  static void								\
  mrc_vec_##type##_create(struct mrc_vec *vec)				\
  {									\
    vec->size_of_type = sizeof(type);					\
  }									\
									\
  static void								\
  mrc_vec_##type##_axpy(struct mrc_vec *y, double alpha, struct mrc_vec *x) \
  {									\
    assert(y->len == x->len);						\
    assert(strcmp(mrc_vec_type(y), mrc_vec_type(x)) == 0);		\
    type *y_arr = y->arr, *x_arr =  x->arr;				\
    type talpha = (type) alpha;						\
    for (int i = 0; i < y->len; i++) {					\
      y_arr[i] += talpha * x_arr[i];					\
     }									\
  }									\
									\
  static void								\
  mrc_vec_##type##_waxpy(struct mrc_vec *w, double alpha, struct mrc_vec *x, struct mrc_vec *y)	\
  {									\
    assert(y->len == x->len);						\
    assert(w->len == x->len);						\
    assert(strcmp(mrc_vec_type(y), mrc_vec_type(x)) == 0);		\
    assert(strcmp(mrc_vec_type(w), mrc_vec_type(x)) == 0);		\
    type *y_arr = y->arr, *x_arr = x->arr, *w_arr =  w->arr;		\
    type talpha = (type) alpha;						\
    for (int i = 0; i < y->len; i++) {					\
      w_arr[i] = talpha * x_arr[i] + y_arr[i];				\
    }									\
  }									\
									\
  static void								\
  mrc_vec_##type##_axpby(struct mrc_vec *y, double alpha, struct mrc_vec *x, double beta) \
  {									\
    assert(y->len == x->len);						\
    assert(strcmp(mrc_vec_type(y), mrc_vec_type(x)) == 0);		\
    type *y_arr = y->arr, *x_arr =  x->arr;				\
    type talpha = (type) alpha, tbeta = (type) beta;			\
    for (int i = 0; i < y->len; i++) {					\
      y_arr[i] = talpha * x_arr[i] + tbeta * y_arr[i];			\
     }									\
  }									\
									\
  static void								\
  mrc_vec_##type##_set(struct mrc_vec *x, double val)			\
  {									\
    type *arr = x->arr;							\
    type tval = (type) val;						\
    for (int i = 0; i < x->len; i++) {					\
      arr[i] = tval;							\
    }									\
  }									\
									\
  static void								\
  mrc_vec_##type##_copy(struct mrc_vec *vec_to, struct mrc_vec *vec_from) \
  {									\
    assert(vec_to->len == vec_from->len);				\
    assert(strcmp(mrc_vec_type(vec_to), mrc_vec_type(vec_from)) == 0);	\
    memcpy(vec_to->arr, vec_from->arr, vec_to->len * sizeof(type));	\
  }									\
  									\
  static struct mrc_obj_method mrc_vec_##type##_methods[] = {		\
    MRC_OBJ_METHOD("axpy", mrc_vec_##type##_axpy),			\
    MRC_OBJ_METHOD("waxpy", mrc_vec_##type##_waxpy),			\
    MRC_OBJ_METHOD("set", mrc_vec_##type##_set),			\
    MRC_OBJ_METHOD("copy", mrc_vec_##type##_copy),			\
    {},									\
  };									\
									\
  static struct mrc_vec_ops mrc_vec_##type##_ops = {			\
    .name                  = #type,					\
    .methods               = mrc_vec_##type##_methods,			\
    .create                = mrc_vec_##type##_create,			\
    .set_array             = mrc_vec_sub_set_array,			\
    .replace_array         = mrc_vec_sub_replace_array,			\
    .get_array             = mrc_vec_sub_get_array,			\
    .put_array             = mrc_vec_sub_put_array,			\
    .axpy                  = mrc_vec_##type##_axpy,			\
    .waxpy                 = mrc_vec_##type##_waxpy,			\
    .axpby                 = mrc_vec_##type##_axpby,			\
    .set                   = mrc_vec_##type##_set,			\
    .copy                  = mrc_vec_##type##_copy,			\
  };

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
#ifdef HAVE_PETSC
  mrc_class_register_subclass(&mrc_class_mrc_vec, &mrc_vec_petsc_ops);
#endif
  
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

