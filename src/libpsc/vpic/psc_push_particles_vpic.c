
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

  // needs E, B
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", EX, HX + 6);
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "vpic", MP_DONT_COPY);

  struct vpic_push_particles *vpushp = psc_push_particles_vpic(push)->vpushp;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;

  vpic_push_particles_prep(vpushp, vmprts, vmflds);

  psc_mfields_put_as(mflds, mflds_base, 0, 0);
  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_push_mprts

static void
psc_push_particles_vpic_push_mprts(struct psc_push_particles *push,
				   struct psc_mparticles *mprts_base,
				   struct psc_mfields *mflds_base)
{
  // needs E, B (not really, because they're already in interpolator), rhob?
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", EX, HX + 6);
  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "vpic", 0);

  struct vpic_push_particles *vpushp = psc_push_particles_vpic(push)->vpushp;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;

  vpic_push_particles_push_mprts(vpushp, vmprts, vmflds);

  // update jf FIXME: rhob too, probably, depending on b.c.
  psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
  psc_mparticles_put_as(mprts, mprts_base, 0);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_stagger_mprts

static void
psc_push_particles_vpic_stagger_mprts(struct psc_push_particles *push,
				      struct psc_mparticles *mprts_base,
				      struct psc_mfields *mflds_base)
{
  // needs E, B
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "vpic", EX, HX + 6);

  struct vpic_push_particles *vpushp = psc_push_particles_vpic(push)->vpushp;
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts_base)->vmprts;

  vpic_push_particles_stagger_mprts(vpushp, vmprts, vmflds);

  // updates no fields
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
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

