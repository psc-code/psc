
#ifndef PSC_FIELDS_AS_C_H
#define PSC_FIELDS_AS_C_H

#include "psc_fields_c.h"

typedef fields_c_real_t fields_real_t;
typedef fields_c_t      fields_t;
#define fields_t_from_psc_fields fields_c_t_from_psc_fields

#define MPI_FIELDS_REAL MPI_FIELDS_C_REAL

#define F3(pf, fldnr, jx,jy,jz) F3_C(pf, fldnr, jx,jy,jz)
#define _F3(fld, m, i,j,k) _F3_C(fld, m, i,j,k)

#define fields_t_set_nan              double_set_nan
#define FIELDS_TYPE                   "c"

#define PSC_FIELDS_AS_C 1

#endif
