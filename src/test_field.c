
#include "psc_testing.h"
#include "psc_push_fields.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>

static void
setup_fields(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(&psc, p, jx, jy, jz, &ix, &iy, &iz);
      f_real xx = 2.*M_PI * ix / psc.domain.gdims[0];
      f_real zz = 2.*M_PI * iz / psc.domain.gdims[2];
      F3_BASE(pf, JXI, jx,jy,jz) = cos(xx) * sin(zz);
      F3_BASE(pf, JYI, jx,jy,jz) = sin(xx) * sin(zz);
      F3_BASE(pf, JZI, jx,jy,jz) = cos(xx) * cos(zz);
    } foreach_3d_g_end;
  }
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct psc_mod_config conf_fortran = {
    .mod_field = "fortran",
  };
  struct psc_mod_config conf_c = {
    .mod_field = "c",
  };

  // test push_field_a

  struct psc_case *_case = psc_create_test_xz(&conf_fortran);
  psc_case_setup(_case);
  mfields_base_t *flds = &psc.flds;
  setup_fields(flds);
  // psc_dump_field(EX, "ex0");
  psc_push_fields_step_a(psc.push_fields, flds);
  // psc_dump_field(EX, "ex1");
  psc_save_fields_ref(&psc, flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xz(&conf_c);
  psc_case_setup(_case);
  setup_fields(flds);
  psc_push_fields_step_a(psc.push_fields, flds);
  // psc_dump_field(EX, "ex2");
  // psc_dump_field(EY, "ey2");
  // psc_dump_field(EZ, "ez2");
  // psc_dump_field(HX, "hx2");
  // psc_dump_field(HY, "hy2");
  // psc_dump_field(HZ, "hz2");
  psc_check_fields_ref(&psc, flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

  // test push_field_a

  _case = psc_create_test_xz(&conf_fortran);
  psc_case_setup(_case);
  setup_fields(flds);
  psc_dump_field(flds, EX, "ex0");
  psc_push_fields_step_b(psc.push_fields, flds);
  psc_dump_field(flds, EX, "ex1");
  psc_save_fields_ref(&psc, flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xz(&conf_c);
  psc_case_setup(_case);
  setup_fields(flds);
  psc_push_fields_step_b(psc.push_fields, flds);
  psc_dump_field(flds, EX, "ex2");
  psc_dump_field(flds, EY, "ey2");
  psc_dump_field(flds, EZ, "ez2");
  psc_dump_field(flds, HX, "hx2");
  psc_dump_field(flds, HY, "hy2");
  psc_dump_field(flds, HZ, "hz2");
  psc_check_fields_ref(&psc, flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

  prof_print();

  MPI_Finalize();
}
