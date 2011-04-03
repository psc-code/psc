
#include "psc_testing.h"
#include "psc_bnd.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>

static void
setup_jx(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(&psc, p, jx, jy, jz, &ix, &iy, &iz);
      f_real xx = 2.*M_PI * ix / psc.domain.gdims[0];
      f_real zz = 2.*M_PI * iz / psc.domain.gdims[2];
      F3_BASE(pf, JXI, jx,jy,jz) = cos(xx) * sin(zz);
    } foreach_3d_g_end;
  }
}

static void
setup_jx_noghost(mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(&psc, p, jx, jy, jz, &ix, &iy, &iz);
      f_real xx = 2.*M_PI * ix / psc.domain.gdims[0];
      f_real zz = 2.*M_PI * iz / psc.domain.gdims[2];
      F3_BASE(pf, JXI, jx,jy,jz) = cos(xx) * sin(zz);
    } foreach_3d_end;
  }
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

  struct psc_case *_case = psc_create_test_xz(&conf_fortran);
  psc_case_setup(_case);
  mfields_base_t *flds = &psc.flds;
  setup_jx(flds);
  //  psc_dump_field(JXI, "jx0");
  psc_bnd_add_ghosts(psc.bnd, flds, JXI, JXI + 1);
  //  psc_dump_field(JXI, "jx1");
  psc_save_fields_ref(&psc, flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xz(&conf_c);
  psc_case_setup(_case);
  setup_jx(flds);
  psc_bnd_add_ghosts(psc.bnd, flds, JXI, JXI + 1);
  //  psc_dump_field(JXI, "jx2");
  psc_check_currents_ref_noghost(&psc, flds, 1e-10);
  psc_case_destroy(_case);

  // test psc_fill_ghosts()

  _case = psc_create_test_xz(&conf_fortran);
  psc_case_setup(_case);
  setup_jx_noghost(flds);
  psc_dump_field(flds, JXI, "jx0");
  psc_bnd_fill_ghosts(psc.bnd, flds, JXI, JXI + 1);
  psc_dump_field(flds, JXI, "jx1");
  psc_save_fields_ref(&psc, flds);
  psc_case_destroy(_case);

  _case = psc_create_test_xz(&conf_c);
  psc_case_setup(_case);
  setup_jx_noghost(flds);
  psc_bnd_fill_ghosts(psc.bnd, flds, JXI, JXI + 1);
  psc_dump_field(flds, JXI, "jx2");
  psc_check_currents_ref(&psc, flds, 1e-10);
  psc_case_destroy(_case);

  prof_print();

  MPI_Finalize();
}
