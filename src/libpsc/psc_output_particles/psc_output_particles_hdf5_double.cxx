
#include "psc_output_particles_private.h"
#include "psc_particles_as_double.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>

#include <psc_output_particles_hdf5_common.cxx>

// ======================================================================
// psc_output_particles: subclass "hdf5_double"

struct psc_output_particles_ops_hdf5_double : psc_output_particles_ops {
  psc_output_particles_ops_hdf5_double() {
    name                  = "hdf5_double";
    size                  = sizeof(struct psc_output_particles_hdf5);
    param_descr           = psc_output_particles_hdf5_descr;
    create                = psc_output_particles_hdf5_create;
    setup                 = psc_output_particles_hdf5_setup;
    destroy               = psc_output_particles_hdf5_destroy;
    run                   = psc_output_particles_hdf5_run;
  }
} psc_output_particles_hdf5_double_ops;
