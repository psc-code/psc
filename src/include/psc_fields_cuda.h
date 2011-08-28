
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include "psc.h"
#include "cuda_wrap.h"

typedef float fields_cuda_real_t;

#define MPI_FIELDS_CUDA_REAL MPI_FLOAT

typedef struct {
  fields_cuda_real_t *h_flds;
  fields_cuda_real_t *d_flds;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  char **name; //> name for each component
} fields_cuda_t;

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define F3_OFF_CUDA(pf, fldnr, jx,jy,jz)				\
  ((((((fldnr)								\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_CUDA(pf, fldnr, jx,jy,jz)		\
  ((pf)->h_flds[F3_OFF_CUDA(pf, fldnr, jx,jy,jz)])

#else

#define F3_CUDA(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_CUDA(pf, fldnr, jx,jy,jz);			\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &((pf)->h_flds[off]);						\
    }))

#endif

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

//void fields_cuda_alloc(fields_cuda_t *pf);
//void fields_cuda_free(fields_cuda_t *pf);
void fields_cuda_get(fields_cuda_t *pf, int mb, int me);
void fields_cuda_put(fields_cuda_t *pf, int mb, int me);
//void fields_cuda_zero(fields_cuda_t *pf, int m);

EXTERN_C void __fields_cuda_get(fields_cuda_t *pf);
EXTERN_C void __fields_cuda_put(fields_cuda_t *pf);

#endif
