
#ifndef PSC_FIELDS_AS_C_H
#define PSC_FIELDS_AS_C_H

#include "psc_fields_c.h"

typedef fields_c_t fields_t;
typedef mfields_c_t mfields_t;

#define F3(fldnr, jx,jy,jz) F3_C(pf, fldnr, jx,jy,jz)

#define psc_mfields_get          psc_mfields_c_get
#define psc_mfields_put  	 psc_mfields_c_put
#define fields_zero         fields_c_zero

#endif
