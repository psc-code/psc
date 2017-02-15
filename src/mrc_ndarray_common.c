
#include "mrc_ndarray.h"

// ----------------------------------------------------------------------
// mrc_ndarray_sub_create

static void
mrc_ndarray_sub_create(struct mrc_ndarray *nd)
{
  nd->data_type = MRC_NT_TYPE;
  nd->size_of_type = sizeof(TYPE);
}

// ----------------------------------------------------------------------
// mrc_ndarray_sub_set

static void
mrc_ndarray_sub_set(struct mrc_ndarray *nd, double _val)
{
  TYPE val = (TYPE) _val;

  // OPT, there'd be a point in choosing an iterator that proceeds
  // according to the underlying layout (right now, the standard iterator
  // will do fortran order)
  // also, if the underlying data is contiguous, we could dispatch to vec,
  // or use a faster iterator
  struct mrc_ndarray_it it;

  mrc_ndarray_it_all(&it, nd);
  for (; !mrc_ndarray_it_done(&it); mrc_ndarray_it_next(&it)) {
    IT_TYPE(&it, TYPE) = val;
  }
}

// ----------------------------------------------------------------------
// mrc_ndarray_sub_copy

static void
mrc_ndarray_sub_copy(struct mrc_ndarray *to, struct mrc_ndarray *from)
{
  assert(mrc_ndarray_data_type(to) == mrc_ndarray_data_type(from));
  assert(mrc_ndarray_same_shape(to, from));
  // FIXME, optimize if both arrays are contiguous
  struct mrc_ndarray_it it_to, it_from;

  mrc_ndarray_it_all(&it_to, to);
  mrc_ndarray_it_all(&it_from, from);
  for (; !mrc_ndarray_it_done(&it_to); mrc_ndarray_it_next(&it_to), mrc_ndarray_it_next(&it_from)) {
    IT_TYPE(&it_to, TYPE) = IT_TYPE(&it_from, TYPE);
  }
}

// ----------------------------------------------------------------------
// mrc_ndarray_sub_ops

struct mrc_ndarray_ops mrc_ndarray_sub_ops = {
  .name                  = mrc_ndarray_sub_name,
  .create                = mrc_ndarray_sub_create,
  .set                   = mrc_ndarray_sub_set,
  .copy                  = mrc_ndarray_sub_copy,
};

