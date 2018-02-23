
#include "psc_push_fields_private.h"

#include "psc_fields_vpic.h"
#include "psc_method.h"

#include "push_fields.hxx"

#include "vpic_iface.h"

// ----------------------------------------------------------------------

class psc_push_fields_vpic : PushFieldsBase
{
public:
  psc_push_fields_vpic()
  {
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim_);
  }

  void push_E(struct psc_push_fields *push, struct psc_mfields *mflds_base,
	      double dt_fac) override
  {
    // needs J, E, B, TCA, material
    mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(JXI, VPIC_MFIELDS_N_COMP);
    FieldArray *vmflds = psc_mfields_vpic(mf.mflds())->vmflds_fields;
    Simulation_push_mflds_E(sim_, vmflds, dt_fac);
    Simulation_field_injection(sim_); // FIXME, this isn't the place, should have its own psc_field_injection
    
    // updates E, TCA, and B ghost points FIXME 9 == TCAX
    mf.put_as(mflds_base, EX, 9 + 3);
  }

  void push_H(struct psc_push_fields *push, struct psc_mfields *mflds_base,
	      double dt_fac) override
  {
    // needs E, B
    mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(EX, HX + 6);
    FieldArray *vmflds = psc_mfields_vpic(mf.mflds())->vmflds_fields;
    Simulation_push_mflds_H(sim_, vmflds, dt_fac);
    
    // updates B
    mf.put_as(mflds_base, HX, HX + 3);
  }

private:
  Simulation *sim_;
};

// ----------------------------------------------------------------------
// psc_push_fields_vpic_setup

static void
psc_push_fields_vpic_setup(struct psc_push_fields *push)
{
  PscPushFields<psc_push_fields_vpic> pushf(push);
  new(pushf.sub()) psc_push_fields_vpic;
}

// ----------------------------------------------------------------------
// psc_push_fields_vpic_destroy

static void
psc_push_fields_vpic_destroy(struct psc_push_fields *push)
{
  PscPushFields<psc_push_fields_vpic> pushf(push);
  pushf.sub()->~psc_push_fields_vpic();
}

// ----------------------------------------------------------------------
// psc_push_fields_vpic_push_mflds_H

static void
psc_push_fields_vpic_push_mflds_H(struct psc_push_fields *push,
				  struct psc_mfields *mflds_base,
				  double dt_fac)
{
  PscPushFields<psc_push_fields_vpic> pushf(push);
  pushf->push_H(push, mflds_base, dt_fac);
}

// ----------------------------------------------------------------------
// psc_push_fields_vpic_push_mflds_E

static void
psc_push_fields_vpic_push_mflds_E(struct psc_push_fields *push,
				  struct psc_mfields *mflds_base,
				  double dt_fac)
{
  PscPushFields<psc_push_fields_vpic> pushf(push);
  pushf->push_E(push, mflds_base, dt_fac);
}

// ----------------------------------------------------------------------
// psc_push_fields: subclass "vpic"

struct psc_push_fields_ops_vpic : psc_push_fields_ops {
  psc_push_fields_ops_vpic() {
    name                  = "vpic";
    size                  = sizeof(struct psc_push_fields_vpic);
    setup                 = psc_push_fields_vpic_setup;
    push_mflds_H          = psc_push_fields_vpic_push_mflds_H;
    push_mflds_E          = psc_push_fields_vpic_push_mflds_E;
  }
} psc_push_fields_vpic_ops;

