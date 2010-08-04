
#include "psc_testing.h"
#include "util/profile.h"
#include "util/params.h"

#include <stdio.h>
#include <math.h>
#include <mpi.h>

static void
setup_jx()
{
  for (int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++) {
    for (int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++) {
      for (int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++) {
	f_real xx = 2.*M_PI * ix / psc.domain.itot[0];
	f_real zz = 2.*M_PI * iz / psc.domain.itot[2];
	F3_BASE(JXI, ix,iy,iz) = cos(xx) * sin(zz);
      }
    }
  }
}

static void
setup_jx_noghost()
{
  for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
    for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
      for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {
	f_real xx = 2.*M_PI * ix / psc.domain.itot[0];
	f_real zz = 2.*M_PI * iz / psc.domain.itot[2];
	F3_BASE(JXI, ix,iy,iz) = cos(xx) * sin(zz);
      }
    }
  }
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  // test psc_add_ghosts()

  struct psc_mod_config conf_fortran = {
    .mod_bnd = "fortran",
  };
  struct psc_mod_config conf_c = {
    .mod_bnd = "c",
  };

  psc_create_test_xz(&conf_fortran);
  setup_jx();
  //  psc_dump_field(JXI, "jx0");
  psc_add_ghosts(JXI, JXI + 1);
  //  psc_dump_field(JXI, "jx1");
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_c);
  setup_jx();
  psc_add_ghosts(JXI, JXI + 1);
  //  psc_dump_field(JXI, "jx2");
  psc_check_currents_ref_noghost(1e-10);
  psc_destroy();

  // test psc_fill_ghosts()

  psc_create_test_xz(&conf_fortran);
  setup_jx_noghost();
  psc_dump_field(JXI, "jx0");
  psc_fill_ghosts(JXI, JXI + 1);
  psc_dump_field(JXI, "jx1");
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_c);
  setup_jx_noghost();
  psc_fill_ghosts(JXI, JXI + 1);
  psc_dump_field(JXI, "jx2");
  psc_check_currents_ref(1e-10);
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
