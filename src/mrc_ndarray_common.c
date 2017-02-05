
#include "mrc_ndarray.h"

static void
mrc_ndarray_sub_create(struct mrc_ndarray *nd)
{
  nd->data_type = MRC_NT_TYPE;
  nd->size_of_type = sizeof(TYPE);
}

struct mrc_ndarray_ops mrc_ndarray_sub_ops = {
  .name                  = mrc_ndarray_sub_name,
  .create                = mrc_ndarray_sub_create,
};

