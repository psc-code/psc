
#include "psc_push_particles_private.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_push_particles_vpic_prep

static void
psc_push_particles_vpic_prep(struct psc_push_particles *push,
			     struct psc_mparticles *mprts_base,
			     struct psc_mfields *mflds_base)
{
  vpic_load_interpolator_array();
}

// ----------------------------------------------------------------------
// psc_push_particles_vpic_run

static void
psc_push_particles_vpic_push_mprts(struct psc_push_particles *push,
				   struct psc_mparticles *mprts_base,
				   struct psc_mfields *mflds_base)
{
  // FIXME, this is kinda too much stuff all in here,
  // so it should be split up, but it'll do for now
  vpic_clear_accumulator_array();
  vpic_advance_p();
  vpic_emitter();
  vpic_reduce_accumulator_array();
  vpic_boundary_p();
  vpic_calc_jf();
  vpic_current_injection();
}

// ----------------------------------------------------------------------
// psc_push_particles: subclass "vpic"

struct psc_push_particles_ops psc_push_particles_vpic_ops = {
  .name                  = "vpic",
  .prep                  = psc_push_particles_vpic_prep,
  .push_mprts            = psc_push_particles_vpic_push_mprts,
};

