
#include "mrc_fld.h"

#include <mrc_vec.h>

// Don't like dirtying main libmrc code in this way
#ifdef HAVE_PETSC
#include <petscconf.h>
#endif

// ======================================================================
// mrc_ndarray

// ----------------------------------------------------------------------
// _mrc_ndarray_destroy

static void
_mrc_ndarray_destroy(struct mrc_ndarray *nd)
{
  if (nd->arr) {
    mrc_vec_put_array(nd->vec, nd->arr);
    nd->arr = NULL;
  }
}

// ----------------------------------------------------------------------
// get_vec_type
//
// helper to get type for our vec member

static const char *
get_vec_type(struct mrc_ndarray *nd)
{
  const char *vec_type;
  switch (nd->data_type) {
  case MRC_NT_FLOAT:  vec_type = "float";  break;
  case MRC_NT_DOUBLE: vec_type = "double"; break;
  case MRC_NT_INT:    vec_type = "int";    break;
  default: assert(0);
  }

  // This dispatch is still ugly, and maybe we should just, e.g., have
  // mrc_vec "float"/"double" use the petsc implementation as appropriate
#if defined(PETSC_USE_REAL_SINGLE) && !defined(PETSC_USE_COMPLEX)
  if (nd->data_type == MRN_NT_FLOAT) {
    vec_type = "petsc";
  }
#endif
#if defined(PETSC_USE_REAL_DOUBLE) && !defined(PETSC_USE_COMPLEX)
  if (nd->data_type == MRN_NT_DOUBLE) {
    vec_type = "petsc";
  }
#endif

  return vec_type;

  // FIXME? we don't have this info here. The call could be done from mrc_fld,
  // but I'm not so sure it's really needed in the first place.
  /* if (fld->_aos && (strcmp(vec_type, "petsc")==0)){ */
  /*   mrc_vec_set_param_int(fld->_nd->vec, "block_size", fld->_nr_comps); */
  /* }   */
}

// ----------------------------------------------------------------------
// _mrc_ndarray_setup

static void
_mrc_ndarray_setup(struct mrc_ndarray *nd)
{
  switch (nd->data_type) {
  case MRC_NT_FLOAT:  nd->size_of_type = sizeof(float);  break;
  case MRC_NT_DOUBLE: nd->size_of_type = sizeof(double); break;
  case MRC_NT_INT:    nd->size_of_type = sizeof(int);    break;
  default: assert(0);
  };
  
  int n_dims = nd->dims.nr_vals;
  assert(n_dims == nd->offs.nr_vals);
  assert(n_dims == nd->perm.nr_vals);
  assert(n_dims <= MRC_FLD_MAXDIMS);
  nd->n_dims = n_dims;

  int *dims = nd->dims.vals, *offs = nd->offs.vals, *perm = nd->perm.vals;

  nd->len = 1;
  for (int d = 0; d < n_dims; d++) {
    nd->len *= dims[d];
  }

  if (!nd->view_base) {
    // new ndarray
    for (int d = 0; d < n_dims; d++) {
      nd->start[d] = offs[d];
      nd->acc.stride[perm[d]] = 1;
      for (int dd = 0; dd < d; dd++) {
	nd->acc.stride[perm[d]] *= dims[perm[dd]];
      }
    }
  } else {
    // make a view
    // FIXME, this doesn't support perm at this point
    struct mrc_ndarray *nd_base = nd->view_base;
    assert(n_dims == nd_base->n_dims);
    assert(n_dims == nd->view_offs.nr_vals);
    int *base_offs = nd_base->offs.vals, *base_dims = nd_base->dims.vals;
    int *view_offs = nd->view_offs.vals;

    for (int d = 0; d < n_dims; d++) {
      assert(view_offs[d] >= offs[d]);
      assert(view_offs[d] + dims[d] <= base_offs[d] + base_dims[d]);
      nd->acc.stride[d] = nd_base->acc.stride[d];
      nd->start[d] = nd_base->start[d] - view_offs[d] + offs[d];
    }
  }

  if (nd->view_base) {
    mrc_vec_set_array(nd->vec, nd->view_base->arr);
  }

  mrc_vec_set_type(nd->vec, get_vec_type(nd));
  mrc_vec_set_param_int(nd->vec, "len", nd->len);
  mrc_vec_setup(nd->vec);

  // set up arr
  nd->arr = mrc_vec_get_array(nd->vec);
  assert(nd->arr);

  // set up arr_off
  int off = 0;
  for (int d = 0; d < MRC_FLD_MAXDIMS; d++) {
    off += nd->start[d] * nd->acc.stride[d];
  }
  nd->acc.arr_off = nd->arr - off * nd->size_of_type;
}

// ----------------------------------------------------------------------
// mrc_ndarray_set_array

void
mrc_ndarray_set_array(struct mrc_ndarray *nd, void *arr)
{
  mrc_vec_set_array(nd->vec, arr);
}

// ----------------------------------------------------------------------
// mrc_ndarray_replace_array

void
mrc_ndarray_replace_array(struct mrc_ndarray *nd, void *arr)
{
  mrc_vec_replace_array(nd->vec, arr);
}

// ----------------------------------------------------------------------
// mrc_ndarray_access
//
// This returns the needed info (struct mrc_ndarray_access) for easy access to
// the data

struct mrc_ndarray_access *
mrc_ndarray_access(struct mrc_ndarray *nd)
{
  return &nd->acc;
}

// ----------------------------------------------------------------------
// mrc_ndarray_init

static struct mrc_ndarray_ops mrc_ndarray_float_ops;
static struct mrc_ndarray_ops mrc_ndarray_double_ops;
static struct mrc_ndarray_ops mrc_ndarray_int_ops;

static void
mrc_ndarray_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_ndarray, &mrc_ndarray_float_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ndarray, &mrc_ndarray_double_ops);
  mrc_class_register_subclass(&mrc_class_mrc_ndarray, &mrc_ndarray_int_ops);
}

// ----------------------------------------------------------------------
// mrc_ndarray class description

#define VAR(x) (void *)offsetof(struct mrc_ndarray, x)
static struct param mrc_ndarray_descr[] = {
  { "offs"            , VAR(offs)           , PARAM_INT_ARRAY(0, 0) },
  { "dims"            , VAR(dims)           , PARAM_INT_ARRAY(0, 0) },
  { "perm"            , VAR(perm)           , PARAM_INT_ARRAY(0, 0) },
  { "view_base"       , VAR(view_base)      , PARAM_OBJ(mrc_ndarray)},
  { "view_offs"       , VAR(view_offs)      , PARAM_INT_ARRAY(0, 0) },

  { "data_type"       , VAR(data_type)      , MRC_VAR_INT           },
  { "size_of_type"    , VAR(size_of_type)   , MRC_VAR_INT           },
  { "len"             , VAR(len)            , MRC_VAR_INT           },
  { "vec"             , VAR(vec)            , MRC_VAR_OBJ(mrc_vec)  },
  {},
};
#undef VAR

struct mrc_class_mrc_ndarray mrc_class_mrc_ndarray = {
  .name         = "mrc_ndarray",
  .size         = sizeof(struct mrc_ndarray),
  .param_descr  = mrc_ndarray_descr,
  .init         = mrc_ndarray_init,
  .setup        = _mrc_ndarray_setup,
  .destroy      = _mrc_ndarray_destroy,
};

// ----------------------------------------------------------------------
// mrc_ndarray: float, double, int subclasses

#define MAKE_MRC_NDARRAY_TYPE(NAME, type, TYPE)				\
									\
  static void								\
  mrc_ndarray_##NAME##_create(struct mrc_ndarray *nd)			\
  {									\
    nd->data_type = MRC_NT_##TYPE;					\
    nd->size_of_type = sizeof(type);					\
  }									\
  									\
  static struct mrc_ndarray_ops mrc_ndarray_##NAME##_ops = {		\
    .name                  = #NAME,					\
    .create                = mrc_ndarray_##NAME##_create,		\
  };									\

MAKE_MRC_NDARRAY_TYPE(float, float, FLOAT)
MAKE_MRC_NDARRAY_TYPE(double, double, DOUBLE)
MAKE_MRC_NDARRAY_TYPE(int, int, INT)

