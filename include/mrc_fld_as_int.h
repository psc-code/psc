
#ifndef MRC_FLD_AS_INT_H
#define MRC_FLD_AS_INT_H

#include <stdlib.h>

typedef int mrc_fld_data_t;
#define mrc_fld_abs abs

static inline int
mrc_fld_min(int a, int b)
{
  return a > b ? b : a;
}

static inline int
mrc_fld_max(int a, int b)
{
  return a > b ? a : b;
}

#define F3(f, m, i,j,k) MRC_I4(f, i,j,k, m)
#define M3(f, m, i,j,k, p) MRC_I5(f, i,j,k, m, p)
#define FLD_TYPE "int"
#define MPI_MRC_FLD_DATA_T MPI_INT

#endif
