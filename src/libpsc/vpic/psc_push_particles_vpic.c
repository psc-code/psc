
#include "psc_push_particles_private.h"

#include "psc_push_particles_vpic.h"
#include "psc_fields_vpic.h"
#include "psc_particles_vpic.h"
#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_push_particles_vpic_setup

static void
psc_push_particles_vpic_setup(struct psc_push_particles *push)
{
  struct psc_push_particles_vpic *sub = psc_push_particles_vpic(push);

  sub->vpushp = vpic_push_particles_create();
  vpic_push_particles_ctor_from_simulation(sub->vpushp);
  psc_push_particles_setup_super(push);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_prep

static void
psc_push_particles_vpic_prep(struct psc_push_particles *push,
			     struct psc_mparticles *mprts_base,
			     struct psc_mfields *mflds_base)
{
  // At end of step:
  // Fields are updated ... load the interpolator for next time step and
  // particle diagnostics in user_diagnostics if there are any particle
  // species to worry about
  struct vpic_push_particles *vpushp = psc_push_particles_vpic(push)->vpushp;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts_base)->vmprts;

  vpic_push_particles_prep(vpushp, vmprts, vmflds);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_push_mprts

static void
psc_push_particles_vpic_push_mprts(struct psc_push_particles *push,
				   struct psc_mparticles *mprts_base,
				   struct psc_mfields *mflds_base)
{
  struct vpic_push_particles *vpushp = psc_push_particles_vpic(push)->vpushp;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts_base)->vmprts;

  vpic_push_particles_push_mprts(vpushp, vmprts, vmflds);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_stagger_mprts

static void
psc_push_particles_vpic_stagger_mprts(struct psc_push_particles *push,
				      struct psc_mparticles *mprts_base,
				      struct psc_mfields *mflds_base)
{
  struct vpic_push_particles *vpushp = psc_push_particles_vpic(push)->vpushp;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds_base)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts_base)->vmprts;

  vpic_push_particles_stagger_mprts(vpushp, vmprts, vmflds);
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "vpic"

struct psc_push_particles_ops psc_push_particles_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_push_particles_vpic),
  .setup                 = psc_push_particles_vpic_setup,
  .prep                  = psc_push_particles_vpic_prep,
  .push_mprts            = psc_push_particles_vpic_push_mprts,
  .stagger_mprts         = psc_push_particles_vpic_stagger_mprts,
};

