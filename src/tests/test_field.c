
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
      psc_local_to_global_indices(&psc, p, jx, jy, jz, &ix, &iy, &iz);
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

  // test push_field_a

  struct psc_case *_case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "fortran");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  // psc_dump_field(EX, "ex0");
  psc_push_fields_step_a(psc.push_fields, psc.flds);
  psc_save_fields_ref(&psc, psc.flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "c");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_a(psc.push_fields, psc.flds);
  psc_check_fields_ref(&psc, psc.flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

#ifdef USE_CBE
  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "cbe");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_a(psc.push_fields, psc.flds);
  psc_check_fields_ref(&psc, psc.flds,
		       (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);
#endif


  // test push_field_b

  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "fortran");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_b(psc.push_fields, psc.flds);
  psc_save_fields_ref(&psc, psc.flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "c");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_b(psc.push_fields, psc.flds);
  psc_check_fields_ref(&psc, psc.flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

  // test push_field_b

  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "fortran");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_b(psc.push_fields, psc.flds);
  psc_save_fields_ref(&psc, psc.flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "c");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_b(psc.push_fields, psc.flds);
  psc_check_fields_ref(&psc, psc.flds,
		       (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

#ifdef USE_CBE
  _case = psc_create_test_xy();
  psc_push_fields_set_type(psc.push_fields, "cbe");
  psc_case_setup(_case);
  setup_fields(psc.flds);
  psc_push_fields_step_b(psc.push_fields, psc.flds);
  psc_check_fields_ref(&psc, psc.flds,
		       (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);
#endif

  prof_print();

  MPI_Finalize();
}
