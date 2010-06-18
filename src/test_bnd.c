
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

static void
dump_jx(const char *pfx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char fname[100];
  sprintf(fname, "%s-p%d.h5", pfx, rank);
  psc_dump_field(JXI, fname);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  struct psc_mod_config conf_fortran = {
    .mod_bnd = "fortran",
  };
  psc_create_test_xz(&conf_fortran);
  setup_jx();
  dump_jx("jx0");
  psc_fax(JXI);
  dump_jx("jx1");
  psc_save_fields_ref();
  psc_destroy();

  struct psc_mod_config conf_c = {
    .mod_bnd = "c",
  };
  psc_create_test_xz(&conf_c);
  setup_jx();
  psc_fax(JXI);
  dump_jx("jx2");
  psc_check_currents_ref_noghost(1e-10);
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
