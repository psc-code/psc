
#include "psc_push_fields_private.h"

#include "psc_fields_vpic.h"
#include "psc_method.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------

struct psc_push_fields_vpic {
  struct Simulation *sim;
};

#define psc_push_fields_vpic(push_fields) mrc_to_subobj(push_fields, struct psc_push_fields_vpic)

// ----------------------------------------------------------------------
// psc_push_fields_vpic_setup

static void
psc_push_fields_vpic_setup(struct psc_push_fields *push_fields)
{
  struct psc_push_fields_vpic *sub = psc_push_fields_vpic(push_fields);
  
  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sub->sim);
}

// ----------------------------------------------------------------------
// psc_push_fields_vpic_push_mflds_H

static void
psc_push_fields_vpic_push_mflds_H(struct psc_push_fields *push,
				  struct psc_mfields *mflds_base,
				  double frac)
{
  // needs E, B
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", EX, HX + 3);

  struct FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;
  vpic_mfields_advance_b(vmflds, frac);

  // updates B
  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields_vpic_push_mflds_E

static void
psc_push_fields_vpic_push_mflds_E(struct psc_push_fields *push,
				  struct psc_mfields *mflds_base,
				  double frac)
{
  struct psc_push_fields_vpic *sub = psc_push_fields_vpic(push);
  
  // needs J, E, B, TCA, material
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", JXI, VPIC_MFIELDS_N_COMP);

  struct FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;
  vpic_mfields_advance_e(vmflds, frac);
  Simulation_field_injection(sub->sim); // FIXME, this isn't the place, should have its own psc_field_injection

  // updates E, TCA, and B ghost points FIXME 9 == TCAX
  psc_mfields_put_as(mflds, mflds_base, EX, 9 + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields: subclass "vpic"

struct psc_push_fields_ops psc_push_fields_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_push_fields_vpic),
  .setup                 = psc_push_fields_vpic_setup,
  .push_mflds_H          = psc_push_fields_vpic_push_mflds_H,
  .push_mflds_E          = psc_push_fields_vpic_push_mflds_E,
};

