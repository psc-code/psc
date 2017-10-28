
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_CUDA
#include "psc_fields_common.h"
#undef FTYPE

// ----------------------------------------------------------------------

struct psc_mfields_cuda {
  struct cuda_mfields *cmflds;
  struct cuda_mfields_bnd *cbnd;
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

#define BND (2)

#define F3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *mflds_prm.mx[2] + ((jz)-mflds_prm.ilg[2]))			\
    *mflds_prm.mx[1] + ((jy)-mflds_prm.ilg[1]))				\
   *mflds_prm.mx[0] + (0-mflds_prm.ilg[0]))

#define F3_DEV_YZ(fldnr,jy,jz) \
  (d_flds)[F3_DEV_OFF_YZ(fldnr, jy,jz)]

#define F3_DEV_OFF_YZ_(fldnr, jy,jz)					\
  ((((fldnr)								\
     *c_mx[2] + ((jz)-c_ilg[2]))					\
    *c_mx[1] + ((jy)-c_ilg[1]))						\
   *(2*BND+1) + (0-(-BND)))

#define F3_DEV_YZ_(fldnr,jy,jz) \
  (d_flds)[F3_DEV_OFF_YZ_(fldnr, jy,jz)]

#endif
