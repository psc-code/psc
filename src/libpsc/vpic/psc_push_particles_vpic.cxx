
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
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", EX, HX + 6);
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  Simulation_push_mprts_prep(sub->sim, vmflds);

  psc_mfields_put_as(mflds, mflds_base, 0, 0);
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
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", EX, HX + 6);
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "vpic", 0);
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;
  Particles *vmprts = psc_mparticles_vpic(mprts)->vmprts;

  Simulation_push_mprts(sub->sim, vmprts, vmflds);

  // update jf FIXME: rhob too, probably, depending on b.c.
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
  psc_mparticles_put_as(mprts, mprts_base, 0);
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

