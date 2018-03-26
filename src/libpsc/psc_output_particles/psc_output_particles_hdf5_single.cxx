
#include "psc_output_particles_private.h"
#include "psc_particles_as_single.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>

#include "psc_output_particles_hdf5_common.cxx"

// ======================================================================
// psc_output_particles: subclass "hdf5_single"

struct psc_output_particles_ops_hdf5_single : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_hdf5>;
  psc_output_particles_ops_hdf5_single() {
    name                  = "hdf5_single";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
    run                   = psc_output_particles_hdf5_run;
  }
} psc_output_particles_hdf5_single_ops;

