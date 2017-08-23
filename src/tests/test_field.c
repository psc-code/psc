
#include "psc_testing.h"
#include "psc_push_fields.h"
#include "psc_fields_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#if 0
// Note: I have changed the test fields to lie in the 
// xy plane instead of xz. For most field implementations
// this shouldn't make any difference. The CBE 2D implementation
// will **only** handle the xy plane, so I need that to be testable.
// --steve 
static void
setup_fields(struct psc_mfields *flds_base)
{
  psc_foreach_patch(ppsc, p) {
    struct psc_fields *pf_base = psc_mfields_get_patch(flds_base, p);
    struct psc_fields *pf = psc_fields_get_as(pf_base, "c", 0, 0);
    psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(ppsc, p, jx, jy, jz, &ix, &iy, &iz);
      f_real xx = 2.*M_PI * ix / ppsc->domain.gdims[0];
      f_real yy = 2.*M_PI * iy / ppsc->domain.gdims[1];
      F3(pf, JXI, jx,jy,jz) = cos(xx) * sin(yy);
      F3(pf, JYI, jx,jy,jz) = sin(xx) * sin(yy);
      F3(pf, JZI, jx,jy,jz) = cos(xx) * cos(yy);
    } foreach_3d_g_end;
    psc_fields_put_as(pf, pf_base, JXI, JXI + 3);
  }
}
#endif

int
main(int argc, char **argv)
{
#if 0
  psc_testing_init(&argc, &argv);

  // test push_field_a

  struct psc_case *_case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "fortran");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  // psc_dump_field(EX, "ex0");
  psc_push_fields_step_a(ppsc->push_fields, ppsc->flds);
  psc_save_fields_ref(ppsc, ppsc->flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "c");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_a(ppsc->push_fields, ppsc->flds);
  psc_check_fields_ref(ppsc, ppsc->flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

#ifdef USE_CBE
  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "cbe");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_a(ppsc->push_fields, ppsc->flds);
  psc_check_fields_ref(ppsc, ppsc->flds,
		       (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);
#endif


  // test push_field_b

  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "fortran");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_b1(ppsc->push_fields, ppsc->flds);
  psc_push_fields_step_b2(ppsc->push_fields, ppsc->flds);
  psc_save_fields_ref(ppsc, ppsc->flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "c");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_b1(ppsc->push_fields, ppsc->flds);
  psc_push_fields_step_b2(ppsc->push_fields, ppsc->flds);
  psc_check_fields_ref(ppsc, ppsc->flds, (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

  // test push_field_b

  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "fortran");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_b1(ppsc->push_fields, ppsc->flds);
  psc_push_fields_step_b2(ppsc->push_fields, ppsc->flds);
  psc_save_fields_ref(ppsc, ppsc->flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "c");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_b1(ppsc->push_fields, ppsc->flds);
  psc_push_fields_step_b2(ppsc->push_fields, ppsc->flds);
  psc_check_fields_ref(ppsc, ppsc->flds,
		       (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);

#ifdef USE_CBE
  _case = psc_create_test_xy();
  psc_push_fields_set_type(ppsc->push_fields, "cbe");
  psc_case_setup(_case);
  setup_fields(ppsc->flds);
  psc_push_fields_step_b1(ppsc->push_fields, ppsc->flds);
  psc_push_fields_step_b2(ppsc->push_fields, ppsc->flds);
  psc_check_fields_ref(ppsc, ppsc->flds,
		       (int []) { EX, EY, EZ, HX, HY, HZ, -1 }, 1e-7);
  psc_case_destroy(_case);
#endif

  psc_testing_finalize();
#endif
}
