
#define TYPE float
#define MRC_NT_TYPE MRC_NT_FLOAT
#define mrc_ndarray_sub_name "float"
#define mrc_ndarray_sub_ops mrc_ndarray_float_ops

// FIXME, this kinda breaks layers as we're just doing mrc_ndarray here,
// but we're only using mrc_fld_max() etc
#include "mrc_fld_as_float.h"

#include "mrc_ndarray_common.c"

