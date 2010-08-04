
#ifndef PSC_FIELDS_AS_C_H
#define PSC_FIELDS_AS_C_H

#include "psc_fields_c.h"

typedef psc_fields_c_t psc_fields_t;

#define F3(fldnr, jx,jy,jz) F3_C(pf, fldnr, jx,jy,jz)

#define psc_fields_get          psc_fields_c_get
#define psc_fields_put  	psc_fields_c_put
#define psc_fields_zero         psc_fields_c_zero

#endif
