
#include <psc.h>

#include <mrc_params.h>

int
psc_main(int *argc, char ***argv, struct psc_ops *type)
{
  MPI_Init(argc, argv);
  libmrc_params_init(*argc, *argv);

  mrc_class_register_subclass(&mrc_class_psc, type);

  // psc_create() will create the psc object, create the sub-objects
  // (particle, field pusher and many others) and set the parameter defaults.
  // It will also set the psc_bubble defaults and call psc_bubble_create(),
  // which will change some of the general defaults to match this case.
  struct psc *psc = psc_create(MPI_COMM_WORLD);

  // psc_set_from_options() will override general and bubble psc parameters
  // if given on the command line. It will also call
  // psc_bubble_set_from_options()
  psc_set_from_options(psc);

  // psc_setup() will set up the various sub-objects (particle pusher, ...)
  // and set up the initial domain partition, the particles and the fields.
  // The standard implementation, used here, will set particles using
  // psc_bubble_init_npt and the fields using setup_field()
  psc_setup(psc);

  // psc_view() will just print a whole lot of info about the psc object and
  // sub-objects, in particular all the parameters.
  psc_view(psc);

  // psc_integrate() uses the standard implementation, which does the regular
  // classic PIC time integration loop
  psc_integrate(psc);

  // psc_destroy() just cleans everything up when we're done.
  psc_destroy(psc);

  libmrc_params_finalize();
  MPI_Finalize();

  return 0;
}
