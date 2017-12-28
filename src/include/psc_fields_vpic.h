
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "psc_fields_private.h"

#include <mrc_common.h>

#define FTYPE FTYPE_VPIC
#include "psc_fields_common.h"
#undef FTYPE

#include "vpic_iface.h"

BEGIN_C_DECLS

// ----------------------------------------------------------------------

struct psc_mfields_vpic {
  Simulation *sim;
  FieldArray *vmflds_fields;
  HydroArray *vmflds_hydro;
};

#define psc_mfields_vpic(mflds) ({					\
      assert((struct psc_mfields_ops *) mflds->obj.ops == &psc_mfields_vpic_ops); \
      mrc_to_subobj(mflds, struct psc_mfields_vpic);			\
    })


// These should otherwise be static, but they are defined in a separate C++ source file

void psc_mfields_vpic_setup(struct psc_mfields *mflds);
void psc_mfields_vpic_destroy(struct psc_mfields *mflds);
fields_vpic_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p);

END_C_DECLS

#endif
