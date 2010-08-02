
#include "psc.h"
#include "util/profile.h"
#include "util/params.h"

#include <stdio.h>
#include <math.h>
#include <mpi.h>

static void
setup_fields()
{
  for (int iz = psc.ilg[2]; iz < psc.ihg[2]; iz++) {
    for (int iy = psc.ilg[1]; iy < psc.ihg[1]; iy++) {
      for (int ix = psc.ilg[0]; ix < psc.ihg[0]; ix++) {
	f_real xx = 2.*M_PI * ix / psc.domain.itot[0];
	f_real zz = 2.*M_PI * iz / psc.domain.itot[2];
	F3_BASE(JXI, ix,iy,iz) = cos(xx) * sin(zz);
	F3_BASE(JYI, ix,iy,iz) = sin(xx) * sin(zz);
	F3_BASE(JZI, ix,iy,iz) = cos(xx) * cos(zz);
      }
    }
  }
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  struct psc_mod_config conf_fortran = {
    .mod_field = "fortran",
  };
  struct psc_mod_config conf_c = {
    .mod_field = "c",
  };

  // test push_field_a

  psc_create_test_xz(&conf_fortran);
  setup_fields();
  // psc_dump_field(EX, "ex0");
  psc_push_field_a();
  // psc_dump_field(EX, "ex1");
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_c);
  setup_fields();
  psc_push_field_a();
  // psc_dump_field(EX, "ex2");
  // psc_dump_field(EY, "ey2");
  // psc_dump_field(EZ, "ez2");
  // psc_dump_field(BX, "bx2");
  // psc_dump_field(BY, "by2");
  // psc_dump_field(BZ, "bz2");
  psc_check_fields_ref((int []) { EX, EY, EZ, BX, BY, BZ, -1 }, 1e-7);
  psc_destroy();

  // test push_field_a

  psc_create_test_xz(&conf_fortran);
  setup_fields();
  psc_dump_field(EX, "ex0");
  psc_push_field_b();
  psc_dump_field(EX, "ex1");
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_c);
  setup_fields();
  psc_push_field_b();
  psc_dump_field(EX, "ex2");
  psc_dump_field(EY, "ey2");
  psc_dump_field(EZ, "ez2");
  psc_dump_field(BX, "bx2");
  psc_dump_field(BY, "by2");
  psc_dump_field(BZ, "bz2");
  psc_check_fields_ref((int []) { EX, EY, EZ, BX, BY, BZ, -1 }, 1e-7);
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
