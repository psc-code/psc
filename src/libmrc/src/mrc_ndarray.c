
#include "mrc_ndarray.h"

#include <mrc_vec.h>
#include <mrc_io.h>

#include <stdio.h>
#include <assert.h>

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

  // offs default is all zeros
  if (nd->offs.nr_vals == 0) {
    mrc_ndarray_set_param_int_array(nd, "offs", n_dims, NULL);
  }

  // perm default is identity map
  if (nd->perm.nr_vals == 0) {
    assert(MRC_NDARRAY_MAXDIMS == 5);
    mrc_ndarray_set_param_int_array(nd, "perm", n_dims, (int []) { 0, 1, 2, 3, 4 });
  }

  assert(n_dims == nd->offs.nr_vals);
  assert(n_dims == nd->perm.nr_vals);
  assert(n_dims <= MRC_NDARRAY_MAXDIMS);
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
      nd->nd_acc.stride[perm[d]] = 1;
      for (int dd = 0; dd < d; dd++) {
	nd->nd_acc.stride[perm[d]] *= dims[perm[dd]];
      }
    }
  } else {
    // make a view
    struct mrc_ndarray *nd_base = nd->view_base;

    // FIXME, this doesn't support perm at this point
    for (int d = 0; d < n_dims; d++) {
      assert(nd->perm.vals[d] == d);
    }

    // view_offs default is offs
    if (nd->view_offs.nr_vals == 0) {
      mrc_ndarray_set_param_int_array(nd, "view_offs", n_dims, offs);
    }

    assert(n_dims == nd_base->n_dims);
    assert(n_dims == nd->view_offs.nr_vals);
    int *base_offs = nd_base->offs.vals, *base_dims = nd_base->dims.vals;
    int *view_offs = nd->view_offs.vals;

    for (int d = 0; d < n_dims; d++) {
      assert(view_offs[d] >= base_offs[d]);
      assert(view_offs[d] + dims[d] <= base_offs[d] + base_dims[d]);
      nd->nd_acc.stride[d] = nd_base->nd_acc.stride[d];
      nd->start[d] = nd_base->start[d] - view_offs[d] + offs[d];
    }
  }

  if (nd->view_base) {
    mrc_vec_set_array(nd->vec, nd->view_base->arr);
  }

  mrc_vec_set_type(nd->vec, get_vec_type(nd));
  assert(nd->len < 0x80000000); // FIXME, should support 64-bit for real
  mrc_vec_set_param_int(nd->vec, "len", nd->len);
  mrc_vec_setup(nd->vec);

  // set up arr
  nd->arr = mrc_vec_get_array(nd->vec);
  assert(nd->arr);

  // set up arr_off
  size_t off = 0;
  for (int d = 0; d < n_dims; d++) {
    off += nd->start[d] * nd->nd_acc.stride[d];
  }
  //mprintf("off %ld\n", off);
  nd->nd_acc.arr_off = nd->arr - off * nd->size_of_type;

  // store more info in nd_acc so we can do bounds checking
  for (int d = 0; d < n_dims; d++) {
    nd->nd_acc.beg[d] = offs[d];
    nd->nd_acc.end[d] = offs[d] + dims[d];
  }
  for (int d = n_dims; d < MRC_NDARRAY_MAXDIMS; d++) {
    nd->nd_acc.beg[d] = 0;
    nd->nd_acc.end[d] = 1;
  }
  nd->nd_acc.data_type = nd->data_type;
}

// ----------------------------------------------------------------------
// mrc_ndarray_read

static void
_mrc_ndarray_read(struct mrc_ndarray *nd, struct mrc_io *io)
{
  // instead of reading back .vec, which wouldn't read its data, anyway,
  // just create a new .vec, then do the usual setup

  nd->vec = mrc_vec_create(mrc_ndarray_comm(nd));
  mrc_ndarray_setup(nd);

  // Now we can actually read the data
  mrc_io_read_ndarray(io, mrc_io_obj_path(io, nd), nd);
}

// ----------------------------------------------------------------------
// mrc_ndarray_write

static void
_mrc_ndarray_write(struct mrc_ndarray *nd, struct mrc_io *io)
{
  mrc_io_write_ndarray(io, mrc_io_obj_path(io, nd), nd);
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
  return &nd->nd_acc;
}

// ----------------------------------------------------------------------
// mrc_ndarray_n_dims

int
mrc_ndarray_n_dims(struct mrc_ndarray *nd)
{
  return nd->n_dims;
}

// ----------------------------------------------------------------------
// mrc_ndarray_dims

int *
mrc_ndarray_dims(struct mrc_ndarray *nd)
{
  return nd->dims.vals;
}

// ----------------------------------------------------------------------
// mrc_ndarray_offs

int *
mrc_ndarray_offs(struct mrc_ndarray *nd)
{
  return nd->offs.vals;
}

// ----------------------------------------------------------------------
// mrc_ndarray_data_type

int
mrc_ndarray_data_type(struct mrc_ndarray *nd)
{
  return nd->data_type;
}

// ----------------------------------------------------------------------
// mrc_ndarray_same_shape
//
// This function checks whether the two ndarrays have compatible shape,
// ie., same number of elements in each dimension, but it doesn't care about
// offsets

bool
mrc_ndarray_same_shape(struct mrc_ndarray *nd1, struct mrc_ndarray *nd2)
{
  if (nd1->n_dims != nd2->n_dims) {
    return false;
  }

  for (int d = 0; d < nd1->n_dims; d++) {
    if (nd1->dims.vals[d] != nd2->dims.vals[d]) {
      return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------
// mrc_ndarray_f_contiguous
//
// checks whether the data layout is contiguous in Fortran order

bool
mrc_ndarray_f_contiguous(struct mrc_ndarray *nd)
{
  size_t stride = 1;
  for (int d = 0; d < nd->n_dims; d++) {
    // if the number of elements in a given dimension is 1, the corresponding
    // stride for that dimension does not matter, so we don't need to check
    if (nd->dims.vals[d] != 1 && nd->nd_acc.stride[d] != stride) {
      return false;
    }
    stride *= nd->dims.vals[d];
  }
  return true;
}

// ----------------------------------------------------------------------
// mrc_ndarray_set

void
mrc_ndarray_set(struct mrc_ndarray *nd, double val)
{
  struct mrc_ndarray_ops *ops = mrc_ndarray_ops(nd);
  assert(ops && ops->set);
  ops->set(nd, val);
}

// ----------------------------------------------------------------------
// mrc_ndarray_copy

void
mrc_ndarray_copy(struct mrc_ndarray *to, struct mrc_ndarray *from)
{
  struct mrc_ndarray_ops *ops = mrc_ndarray_ops(to);
  assert(ops && ops->copy);
  ops->copy(to, from);
}

// ----------------------------------------------------------------------
// mrc_ndarray_scale

void
mrc_ndarray_scale(struct mrc_ndarray *nd, double val)
{
  struct mrc_ndarray_ops *ops = mrc_ndarray_ops(nd);
  assert(ops && ops->scale);
  ops->scale(nd, val);
}

// ----------------------------------------------------------------------
// mrc_ndarray_norm
//
// Note: This returns the local norm, ie., it's not aware of what happens
// on other processors.
// FIXME: should we enforce the mrc_ndarray has size(comm) == 1?

double
mrc_ndarray_norm(struct mrc_ndarray *nd)
{
  struct mrc_ndarray_ops *ops = mrc_ndarray_ops(nd);
  assert(ops && ops->norm);
  return ops->norm(nd);
}

// ----------------------------------------------------------------------
// mrc_ndarray_init

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
  .read         = _mrc_ndarray_read,
  .write        = _mrc_ndarray_write,
  .destroy      = _mrc_ndarray_destroy,
};
