
#include <psc.h>

#ifdef USE_VPIC
#include "vpic/vpic_iface.h"
#endif

#include <mrc_params.h>

int
psc_main(int *argc, char ***argv, struct psc_ops *type)
{
#ifdef USE_VPIC
  vpic_base_init(argc, argv);
#else
  MPI_Init(argc, argv);
#endif
  //  MPI_Init(argc, argv);
  libmrc_params_init(*argc, *argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  mrc_class_register_subclass(&mrc_class_psc, type);

  int from_checkpoint = -1;
  mrc_params_get_option_int("from_checkpoint", &from_checkpoint);

  struct psc *psc;

  if (from_checkpoint < 0) {
    // regular start-up (not from checkpoint)

    // psc_create() will create the psc object, create the sub-objects
    // (particle, field pusher and many others) and set the parameter defaults.
    // It will also set the psc subtype defaults and call psc_subtype_create(),
    // which will change some of the general defaults to match this case.
    psc = psc_create(MPI_COMM_WORLD);
    
    // psc_set_from_options() will override general and bubble psc parameters
    // if given on the command line. It will also call
    // psc_bubble_set_from_options()
    psc_set_from_options(psc);
    
    // psc_setup() will set up the various sub-objects (particle pusher, ...)
    // and set up the initial domain partition, the particles and the fields.
    // The standard implementation, used here, will set particles using
    // psc_bubble_init_npt and the fields using setup_field()
    psc_setup(psc);
  } else {
    // get psc object from checkpoint file

    psc = psc_read_checkpoint(MPI_COMM_WORLD, from_checkpoint);

    int checkpoint_nmax = -1;
    mrc_params_get_option_int("checkpoint_nmax", &checkpoint_nmax);
    if (checkpoint_nmax >= 0) {
      psc->prm.nmax = checkpoint_nmax;
    }

    int checkpoint_wallclock_limit = -1;
    mrc_params_get_option_int("checkpoint_wallclock_limit", &checkpoint_wallclock_limit);
    if (checkpoint_wallclock_limit >= 0) {
      psc->prm.wallclock_limit = checkpoint_wallclock_limit;
    }
  }

  // psc_view() will just print a whole lot of info about the psc object and
  // sub-objects, in particular all the parameters.
  psc_view(psc);
  psc_mparticles_view(psc->particles);
  psc_mfields_view(psc->flds);
  
  // psc_integrate() uses the standard implementation, which does the regular
  // classic PIC time integration loop
  psc_integrate(psc);
  
  // psc_destroy() just cleans everything up when we're done.
  psc_destroy(psc);
  
  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
