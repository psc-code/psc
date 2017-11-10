
#include "psc_push_fields_private.h"

#include "psc_fields_vpic.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_push_fields_vpic_push_mflds_H

static void
psc_push_fields_vpic_push_mflds_H(struct psc_push_fields *push,
				  struct psc_mfields *mflds_base,
				  double frac)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;
  
  vpic_push_fields_advance_b(vmflds, frac);
}

// ----------------------------------------------------------------------
// psc_push_fields_vpic_push_mflds_E

static void
psc_push_fields_vpic_push_mflds_E(struct psc_push_fields *push,
				  struct psc_mfields *mflds_base,
				  double frac)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;

  vpic_push_fields_advance_e(vmflds, frac);
  vpic_field_injection(); // FIXME, this isn't the place, should have its own psc_field_injection
}

// ----------------------------------------------------------------------
// psc_push_fields: subclass "vpic"

struct psc_push_fields_ops psc_push_fields_vpic_ops = {
  .name                  = "vpic",
  .push_mflds_H          = psc_push_fields_vpic_push_mflds_H,
  .push_mflds_E          = psc_push_fields_vpic_push_mflds_E,
};

