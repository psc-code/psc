
#include "psc_testing.h"
#include "psc_push_fields.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>


// Note: I have changed the test fields to lie in the 
// xy plane instead of xz. For most field implementations
// this shouldn't make any difference. The CBE 2D implementation
// will **only** handle the xy plane, so I need that to be testable.
// --steve 
static void
setup_fields(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(p, jx, jy, jz, &ix, &iy, &iz);
      f_real xx = 2.*M_PI * ix / psc.domain.gdims[0];
      f_real yy = 2.*M_PI * iy / psc.domain.gdims[1];
      F3_BASE(pf, JXI, jx,jy,jz) = cos(xx) * sin(yy);
      F3_BASE(pf, JYI, jx,jy,jz) = sin(xx) * sin(yy);
      F3_BASE(pf, JZI, jx,jy,jz) = cos(xx) * cos(yy);
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

#ifdef USE_CBE
  struct psc_mod_config conf_cbe = {
    .mod_field = "cbe",
  };
#endif

  // test push_field_a


  psc_create_test_xy(&conf_fortran);
  mfields_base_t *flds = &psc.flds;
  setup_fields(flds);
  // psc_dump_field(EX, "ex0");
  psc_push_fields_step_a(psc.push_fields, flds);
  // psc_dump_field(EX, "ex1");
  psc_save_fields_ref(flds);
  psc_destroy();

  psc_create_test_xy(&conf_c);
  setup_fields(flds);
  psc_push_fields_step_a(psc.push_fields, flds);
  // psc_dump_field(EX, "ex2");
  // psc_dump_field(EY, "ey2");
  // psc_dump_field(EZ, "ez2");
  // psc_dump_field(HX, "hx2");
  // psc_dump_field(HY, "hy2");
  // psc_dump_field(HZ, "hz2");
  psc_check_fields_ref(flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_destroy();

#ifdef USE_CBE
  psc_create_test_xy(&conf_cbe);
  setup_fields(flds);
  //psc_dump_field(flds,EX,"ex0");
  psc_push_fielda_step_a(psc.push_fields, flds);
  //psc_dump_field(flds,EX, "ex2");
  // psc_dump_field(EY, "ey2");
  // psc_dump_field(EZ, "ez2");
  // psc_dump_field(HX, "hx2");
  // psc_dump_field(HY, "hy2");
  // psc_dump_field(HZ, "hz2");
  psc_check_fields_ref(flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_destroy();
#endif


  // test push_field_b

  psc_create_test_xy(&conf_fortran);
  setup_fields(flds);
  psc_dump_field(flds, EX, "ex0");
  psc_push_fields_step_b(psc.push_fields, flds);
  psc_dump_field(flds, EX, "ex1");
  psc_save_fields_ref(flds);
  psc_destroy();

  psc_create_test_xy(&conf_c);
  setup_fields(flds);
  psc_push_fields_step_b(psc.push_fields, flds);
  psc_dump_field(flds, EX, "ex2");
  psc_dump_field(flds, EY, "ey2");
  psc_dump_field(flds, EZ, "ez2");
  psc_dump_field(flds, HX, "hx2");
  psc_dump_field(flds, HY, "hy2");
  psc_dump_field(flds, HZ, "hz2");
  psc_check_fields_ref(flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_destroy();

#ifdef USE_CBE
  psc_create_test_xy(&conf_cbe);
  setup_fields(flds);
  psc_push_fields_step_b(psc.push_fields, flds);
  //  psc_dump_field(flds,EX, "ex2");
  // psc_dump_field(EY, "ey2");
  // psc_dump_field(EZ, "ez2");
  // psc_dump_field(HX, "hx2");
  // psc_dump_field(HY, "hy2");
  // psc_dump_field(HZ, "hz2");
  psc_check_fields_ref(flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_destroy();
#endif


  prof_print();

  MPI_Finalize();
}
