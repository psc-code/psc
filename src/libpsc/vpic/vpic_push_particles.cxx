
#include "vpic_push_particles.h"

// ======================================================================
// vpic_push_particles implementation

// ----------------------------------------------------------------------
// ctor

vpic_push_particles::vpic_push_particles(Simulation *sim)
  : sim_(sim)
{
}

// ======================================================================
// wrappers

// ----------------------------------------------------------------------
// vpic_push_particles_new_from_Simulation

struct vpic_push_particles *
vpic_push_particles_new_from_Simulation(Simulation *sim)
{
  return new vpic_push_particles(sim);
}

// ----------------------------------------------------------------------
// vpic_push_particles forwards

void vpic_push_particles_push_mprts(struct vpic_push_particles *vpushp,
				    Particles *vmprts, FieldArray *vmflds)
{
  vpushp->sim_->push_mprts(*vmprts, *vmflds);
}

void vpic_push_particles_prep(struct vpic_push_particles *vpushp,
			      Particles *vmprts, FieldArray *vmflds)
{
  vpushp->sim_->push_mprts_prep(*vmflds);
}

