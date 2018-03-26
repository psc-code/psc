
#include "psc_output_particles_private.h"
#include "output_particles.hxx"

struct psc_output_particles_none : PscOutputParticlesParams, OutputParticlesBase
{
  psc_output_particles_none(const PscOutputParticlesParams& params)
    : PscOutputParticlesParams(params)
  {}
};

// ----------------------------------------------------------------------
// psc_output_particles_none_run

static void
psc_output_particles_none_run(struct psc_output_particles *out,
				 struct psc_mparticles *particles_base)
{
}

// ======================================================================
// psc_output_particles: subclass "none"

struct psc_output_particles_ops_none : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_none>;
  psc_output_particles_ops_none() {
    name                  = "none";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
    run                   = psc_output_particles_none_run;
  }
} psc_output_particles_none_ops;
