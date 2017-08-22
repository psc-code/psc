
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_C
#include "psc_fields_common.h"
#undef FTYPE

#define F3_OFF_C(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr - (pf)->first_comp)					\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_C(pf, fldnr, jx,jy,jz)		\
  (((fields_c_real_t *) (pf)->data)[F3_OFF_C(pf, fldnr, jx,jy,jz)])

#else

#define F3_C(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_C(pf, fldnr, jx,jy,jz);				\
      assert(fldnr >= (pf)->first_comp && fldnr < (pf)->first_comp + (pf)->nr_comp); \
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(((fields_c_real_t *) (pf)->data)[off]);				\
    }))

#endif

#endif
