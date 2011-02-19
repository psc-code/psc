
#include "psc_testing.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>

static void
setup_jx()
{
  foreach_3d_g(ix, iy, iz) {
    f_real xx = 2.*M_PI * ix / psc.domain.itot[0];
    f_real zz = 2.*M_PI * iz / psc.domain.itot[2];
    F3_BASE(JXI, ix,iy,iz) = cos(xx) * sin(zz);
  } foreach_3d_g_end;
}

static void
setup_jx_noghost()
{
  foreach_3d(ix, iy, iz, 0, 0) {
    f_real xx = 2.*M_PI * ix / psc.domain.itot[0];
    f_real zz = 2.*M_PI * iz / psc.domain.itot[2];
    F3_BASE(JXI, ix,iy,iz) = cos(xx) * sin(zz);
  } foreach_3d_end;
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

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
  psc_add_ghosts(&psc.pf, JXI, JXI + 1);
  //  psc_dump_field(JXI, "jx1");
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_c);
  setup_jx();
  psc_add_ghosts(&psc.pf, JXI, JXI + 1);
  //  psc_dump_field(JXI, "jx2");
  psc_check_currents_ref_noghost(1e-10);
  psc_destroy();

  // test psc_fill_ghosts()

  psc_create_test_xz(&conf_fortran);
  setup_jx_noghost();
  psc_dump_field(JXI, "jx0");
  psc_fill_ghosts(&psc.pf, JXI, JXI + 1);
  psc_dump_field(JXI, "jx1");
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_c);
  setup_jx_noghost();
  psc_fill_ghosts(&psc.pf, JXI, JXI + 1);
  psc_dump_field(JXI, "jx2");
  psc_check_currents_ref(1e-10);
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
