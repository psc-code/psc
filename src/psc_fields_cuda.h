
#ifndef PSC_FIELD_CUDA_H
#define PSC_FIELD_CUDA_H

#include "psc.h"

typedef float fields_cuda_real_t;

#define MPI_FIELDS_CUDA_REAL MPI_FLOAT

typedef struct {
  fields_cuda_real_t *flds;
  fields_cuda_real_t *d_flds;
} fields_cuda_t;

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define F3_OFF_CUDA(fldnr, jx,jy,jz)					\
  (((((fldnr								\
       *psc.img[2] + ((jz)-psc.ilg[2]))					\
      *psc.img[1] + ((jy)-psc.ilg[1]))					\
     *psc.img[0] + ((jx)-psc.ilg[0]))))

#if 1

#define F3_CUDA(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_CUDA(fldnr, jx,jy,jz)])

#else
//out of range debugging
#define F3_CUDA(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_CUDA(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < NR_FIELDS * psc.fld_size);				\
      &((pf)->flds[off]);						\
    }))

#endif

// ----------------------------------------------------------------------
// macros to access fields from CUDA (device-side)

#define F3_OFF(fldnr, jx,jy,jz)						\
  ((((fldnr)								\
     *d_mx[2] + ((jz)-d_iglo[2]))					\
    *d_mx[1] + ((jy)-d_iglo[1]))					\
   *d_mx[0] + ((jx)-d_iglo[0]))

#if 1

#define F3(fldnr, jx,jy,jz) \
  (d_flds)[F3_OFF(fldnr, jx,jy,jz)]

#else

#define F3(fldnr, jx,jy,jz)						\
  (*({int off = F3_OFF(fldnr, jx,jy,jz);				\
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

#endif
