
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include "psc.h"
#include "cuda_wrap.h"

typedef float fields_cuda_real_t;

#define MPI_FIELDS_CUDA_REAL MPI_FLOAT

typedef struct {
  fields_cuda_real_t *d_flds;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
} fields_cuda_t;

// ----------------------------------------------------------------------
// macros to access fields from CUDA (device-side)

#define F3_DEV_OFF(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *d_mx[2] + ((jz)-d_ilg[2]))					\
    *d_mx[1] + ((jy)-d_ilg[1]))						\
   *d_mx[0] + ((jx)-d_ilg[0]))

#if 1

#define F3_DEV(fldnr, jx,jy,jz) \
  (d_flds)[F3_DEV_OFF(fldnr, jx,jy,jz)]

#else

#define F3_DEV(fldnr, jx,jy,jz)						\
  (*({int off = F3_DEV_OFF(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * d_mx[0] * d_mx[1] * d_mx[2]);		\
      &(d_flds[off]);							\
    }))

#endif

#endif
