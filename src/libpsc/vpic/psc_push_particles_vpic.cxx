
#include "psc_push_particles_private.h"

#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "psc_method.h"
#include "vpic_iface.h"

// ======================================================================
// psc_push_particles_vpic

struct psc_push_particles_vpic {
  Simulation *sim;
};

#define psc_push_particles_vpic(push) mrc_to_subobj(push, struct psc_push_particles_vpic)

// ----------------------------------------------------------------------
// psc_push_particles_vpic_setup

static void
psc_push_particles_vpic_setup(struct psc_push_particles *push)
{
  struct psc_push_particles_vpic *sub = psc_push_particles_vpic(push);

  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sub->sim);

  psc_push_particles_setup_super(push);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_prep

static void
psc_push_particles_vpic_prep(struct psc_push_particles *push,
			     struct psc_mparticles *mprts_base,
			     struct psc_mfields *mflds_base)
{
  struct psc_push_particles_vpic *sub = psc_push_particles_vpic(push);

  // needs E, B
  mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(EX, HX + 6);
  FieldArray *vmflds = psc_mfields_vpic(mf.mflds())->vmflds_fields;

  Simulation_push_mprts_prep(sub->sim, vmflds);

  mf.put_as(mflds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_push_mprts

static void
psc_push_particles_vpic_push_mprts(struct psc_push_particles *push,
				   struct psc_mparticles *mprts_base,
				   struct psc_mfields *mflds_base)
{
  struct psc_push_particles_vpic *sub = psc_push_particles_vpic(push);

  // needs E, B (not really, because they're already in interpolator), rhob?
  mfields_vpic_t mf = mflds_base->get_as<mfields_vpic_t>(EX, HX + 6);
  FieldArray *vmflds = psc_mfields_vpic(mf.mflds())->vmflds_fields;
  mparticles_vpic_t mprts = mprts_base->get_as<mparticles_vpic_t>();
  Particles *vmprts = psc_mparticles_vpic(mprts.mprts())->vmprts;

  Simulation_push_mprts(sub->sim, vmprts, vmflds);

  // update jf FIXME: rhob too, probably, depending on b.c.
  mf.put_as(mflds_base, JXI, JXI + 3);
  mprts.put_as(mprts_base);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "vpic"

struct psc_push_particles_ops psc_push_particles_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_push_particles_vpic),
  .setup                 = psc_push_particles_vpic_setup,
  .prep                  = psc_push_particles_vpic_prep,
  .push_mprts            = psc_push_particles_vpic_push_mprts,
};

