
#define TYPE double
#define MRC_NT_TYPE MRC_NT_DOUBLE
#define mrc_ndarray_sub_name "double"
#define mrc_ndarray_sub_ops mrc_ndarray_double_ops

// FIXME, this kinda breaks layers as we're just doing mrc_ndarray here,
// but we're only using mrc_fld_max() etc
#include "mrc_fld_as_double.h"

#include "mrc_ndarray_common.c"

