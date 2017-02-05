
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
// mrc_ndarray_sub_ops

struct mrc_ndarray_ops mrc_ndarray_sub_ops = {
  .name                  = mrc_ndarray_sub_name,
  .create                = mrc_ndarray_sub_create,
  .set                   = mrc_ndarray_sub_set,
};

