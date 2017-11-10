
#include "psc_push_particles_private.h"

#include "psc_fields_vpic.h"
#include "vpic_iface.h"

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
  vpic_load_interpolator_array(psc_mfields_vpic(mflds_base)->vmflds);
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_run

static void
psc_push_particles_vpic_push_mprts(struct psc_push_particles *push,
				   struct psc_mparticles *mprts_base,
				   struct psc_mfields *mflds_base)
{
  struct psc_mfields_vpic *sub = psc_mfields_vpic(mflds_base);
  
  // FIXME, this is kinda too much stuff all in here,
  // so it should be split up, but it'll do for now
  vpic_clear_accumulator_array();
  vpic_advance_p();
  vpic_emitter();
  vpic_reduce_accumulator_array();
  vpic_boundary_p(sub->vmflds);
  vpic_calc_jf(sub->vmflds);
  vpic_current_injection();
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "vpic"

struct psc_push_particles_ops psc_push_particles_vpic_ops = {
  .name                  = "vpic",
  .prep                  = psc_push_particles_vpic_prep,
  .push_mprts            = psc_push_particles_vpic_push_mprts,
};

