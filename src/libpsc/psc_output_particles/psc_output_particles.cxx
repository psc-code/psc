
#include "psc_output_particles_private.h"

// ======================================================================
// psc_output_particles_init

extern struct psc_output_particles_ops psc_output_particles_hdf5_single_ops;
extern struct psc_output_particles_ops psc_output_particles_hdf5_double_ops;
extern struct psc_output_particles_ops psc_output_particles_none_ops;
extern struct psc_output_particles_ops psc_output_particles_ascii_ops;

static void
psc_output_particles_init()
{
#ifdef HAVE_HDF5
  mrc_class_register_subclass(&mrc_class_psc_output_particles,
                              &psc_output_particles_hdf5_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_particles,
                              &psc_output_particles_hdf5_double_ops);
#endif
  mrc_class_register_subclass(&mrc_class_psc_output_particles,
			      &psc_output_particles_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_particles,
                              &psc_output_particles_ascii_ops);
}

// ======================================================================
// psc_output_particles class

struct mrc_class_psc_output_particles_ : mrc_class_psc_output_particles {
  mrc_class_psc_output_particles_() {
    name             = "psc_output_particles";
    size             = sizeof(struct psc_output_particles);
    init             = psc_output_particles_init;
  }
} mrc_class_psc_output_particles;

