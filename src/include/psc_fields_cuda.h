
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include "psc_fields_private.h"
#include "cuda_wrap.h"

typedef float fields_cuda_real_t;

#define MPI_FIELDS_CUDA_REAL MPI_FLOAT

struct psc_fields_cuda_bnd {
  fields_cuda_real_t *arr_off;
  int im[3];
  int ib[3];
  fields_cuda_real_t *arr;
};

struct psc_fields_cuda {
  fields_cuda_real_t *d_flds;
  fields_cuda_real_t *d_bnd_buf, *h_bnd_buf;
  struct psc_fields_cuda_bnd bnd;
};

#define psc_fields_cuda(pf) mrc_to_subobj(pf, struct psc_fields_cuda)

// ----------------------------------------------------------------------

struct psc_mfields_cuda {
  fields_cuda_real_t *d_flds;
};

#define psc_mfields_cuda(pf) mrc_to_subobj(pf, struct psc_mfields_cuda)

// ----------------------------------------------------------------------
// macros to access fields from CUDA (device-side)

#define F3_DEV_OFF(fldnr, jx,jy,jz)					\
  ((((fldnr)								\
     *d_consts.mx[2] + ((jz)-d_consts.ilg[2]))				\
    *d_consts.mx[1] + ((jy)-d_consts.ilg[1]))				\
   *d_consts.mx[0] + ((jx)-d_consts.ilg[0]))

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

//OPT: precalc offset, don't do ghost points in invar dir

#define F3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *prm.mx[2] + ((jz)-prm.ilg[2]))					\
    *prm.mx[1] + ((jy)-prm.ilg[1]))					\
   *prm.mx[0] + (0-prm.ilg[0]))

#define F3_DEV_YZ(fldnr,jy,jz) \
  (d_flds)[F3_DEV_OFF_YZ(fldnr, jy,jz)]

#endif
