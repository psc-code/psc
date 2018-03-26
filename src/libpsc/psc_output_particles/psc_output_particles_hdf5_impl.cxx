
#include "psc_particles_single.h"
#include "psc_particles_double.h"

#include "output_particles_hdf5_impl.hxx"

// ======================================================================
// psc_output_particles: subclass "hdf5_single"

struct psc_output_particles_ops_hdf5_single : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_hdf5<MparticlesSingle>>;
  psc_output_particles_ops_hdf5_single() {
    name                  = "hdf5_single";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
  }
} psc_output_particles_hdf5_single_ops;

// ======================================================================
// psc_output_particles: subclass "hdf5_double"

struct psc_output_particles_ops_hdf5_double : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_hdf5<MparticlesDouble>>;
  psc_output_particles_ops_hdf5_double() {
    name                  = "hdf5_double";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
  }
} psc_output_particles_hdf5_double_ops;
