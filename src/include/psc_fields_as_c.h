
#ifndef PSC_FIELDS_AS_C_H
#define PSC_FIELDS_AS_C_H

#include "psc_fields_c.h"

typedef fields_c_real_t fields_real_t;
#define MPI_FIELDS_REAL MPI_FIELDS_C_REAL

#define F3(pf, fldnr, jx,jy,jz) F3_C(pf, fldnr, jx,jy,jz)

#define psc_mfields_get_cf            psc_mfields_get_c
#define psc_mfields_put_cf  	      psc_mfields_put_c
#define fields_t_set_nan              double_set_nan
#define FIELDS_TYPE                   "c"

#define PSC_FIELDS_AS_C 1

#endif
