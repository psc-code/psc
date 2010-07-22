
#include "psc.h"
#include "util/profile.h"
#include "util/params.h"

#include <stdio.h>
#include <mpi.h>

static void
setup_jx()
{
  for (int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++) {
    for (int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++) {
      for (int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++) {
	FF3(JXI, ix,iy,iz) = 1.;
      }
    }
  }
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  params_init(argc, argv);

  struct psc_mod_config conf_fortran = {
  };
  psc_create_test_xz(&conf_fortran);
  setup_jx();

  char fname[100];
  sprintf(fname, "jx0-p%d.h5", rank);
  psc_dump_field(JXI, fname);
  psc_fax(JXI);
  sprintf(fname, "jx1-p%d.h5", rank);
  psc_dump_field(JXI, fname);

  psc_save_fields_ref();
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
