
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "psc_fields_private.h"

#include <mrc_common.h>

#define FTYPE FTYPE_VPIC
#include "psc_fields_common.h"
#undef FTYPE

using mfields_vpic_t = mfields3d<fields_vpic_t>;

template<>
inline fields_vpic_t mfields_vpic_t::operator[](int p)
{
  fields_vpic_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p);
  return psc_mfields_vpic_get_field_t(mflds_, p);
}

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


END_C_DECLS

#endif
