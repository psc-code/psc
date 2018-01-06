
#ifndef PSC_FIELDS_ACC_H
#define PSC_FIELDS_ACC_H

#include "psc_fields_private.h"

typedef float fields_acc_real_t;

struct psc_fields_acc {
};

#define psc_fields_acc(flds) mrc_to_subobj(flds, struct psc_fields_acc)

// ----------------------------------------------------------------------
// macros to access C (host) versions of the fields

#define F3_OFF_ACC(flds, fldnr, jx,jy,jz)				\
  ((((((fldnr)								\
       * (flds)->im[2] + ((jz)-(flds)->ib[2]))				\
      * (flds)->im[1] + ((jy)-(flds)->ib[1]))				\
     * (flds)->im[0] + ((jx)-(flds)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_ACC(flds, fldnr, jx,jy,jz)		\
  (((fields_acc_real_t *) (flds)->data)[F3_OFF_ACC(flds, fldnr, jx,jy,jz)])

#else

#define F3_ACC(flds, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_ACC(flds, fldnr, jx,jy,jz);			\
      assert(fldnr >= 0 && fldnr < (flds)->nr_comp);			\
      assert(jx >= (flds)->ib[0] && jx < (flds)->ib[0] + (flds)->im[0]); \
      assert(jy >= (flds)->ib[1] && jy < (flds)->ib[1] + (flds)->im[1]); \
      assert(jz >= (flds)->ib[2] && jz < (flds)->ib[2] + (flds)->im[2]); \
      &(((fields_acc_real_t *) (flds)->data)[off]);			\
    }))

#endif

// ----------------------------------------------------------------------

struct psc_mfields_acc {
  fields_acc_real_t *flds;
  int ib[3], im[3];
};

#define psc_mfields_acc(mflds) mrc_to_subobj(mflds, struct psc_mfields_acc)

#endif
