
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

#define VAR(x) (void *)offsetof(struct psc_output_particles, params.x)
static struct param psc_output_particles_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "basename"           , VAR(basename)             , PARAM_STRING("prt")     },
  { "every_step"         , VAR(every_step)           , PARAM_INT(-1)           },
  { "lo"                 , VAR(lo)                   , PARAM_INT3(0, 0, 0)     },
  { "hi"                 , VAR(hi)                   , PARAM_INT3(0, 0, 0)     },
  { "use_independent_io" , VAR(use_independent_io)   , PARAM_BOOL(false)       },
  { "romio_cb_write"     , VAR(romio_cb_write)       , PARAM_STRING(NULL)      },
  { "romio_ds_write"     , VAR(romio_ds_write)       , PARAM_STRING(NULL)      },
  {},
};
#undef VAR

struct mrc_class_psc_output_particles_ : mrc_class_psc_output_particles {
  mrc_class_psc_output_particles_() {
    name             = "psc_output_particles";
    size             = sizeof(struct psc_output_particles);
    param_descr      = psc_output_particles_descr;
    init             = psc_output_particles_init;
  }
} mrc_class_psc_output_particles;

