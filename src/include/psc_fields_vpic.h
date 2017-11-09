
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_VPIC
#include "psc_fields_common.h"
#undef FTYPE

// ----------------------------------------------------------------------

struct psc_mfields_vpic {
};

#define psc_mfields_vpic(pf) mrc_to_subobj(pf, struct psc_mfields_vpic)

#endif
