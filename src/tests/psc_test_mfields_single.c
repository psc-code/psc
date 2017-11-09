
#include "psc_fields.h"

#include "psc.h"

// ----------------------------------------------------------------------
// test_create_destroy

static void
test_create_destroy(MPI_Comm comm)
{
  struct psc_mfields *mflds = psc_mfields_create(comm);
  psc_mfields_destroy(mflds);
  
}

// ----------------------------------------------------------------------
// dummy_psc_main
//
// FIXME, just a hack

void
dummy_psc_main()
{
  // FIXME: If we don't reference psc_main, we get link errors
  // from missing global variables (for timing etc),
  psc_main(NULL, NULL, NULL);
}

// ----------------------------------------------------------------------
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  test_create_destroy(comm);
  
  MPI_Finalize();
}
