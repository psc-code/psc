
#ifndef PSC_FIELDS_AS_SINGLE_H
#define PSC_FIELDS_AS_SINGLE_H

#include "psc_fields_single.h"

typedef fields_single_real_t fields_real_t;
#define MPI_FIELDS_REAL MPI_FIELDS_SINGLE_REAL

#define F3(pf, fldnr, jx,jy,jz) F3_S(pf, fldnr, jx,jy,jz)

#define fields_t_set_nan              float_set_nan
#define FIELDS_TYPE                   "single"

#define PSC_FIELDS_AS_SINGLE 1

#endif
