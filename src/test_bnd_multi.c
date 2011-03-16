
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
    foreach_3d(p, jx, jy, jz, 0, 0) {
      int ix, iy, iz;
      psc_local_to_global_indices(p, jx, jy, jz, &ix, &iy, &iz);
      F3_BASE(pf, JXI, jx,jy,jz) = iz * 10000 + iy * 100 + ix;
    } foreach_3d_end;
  }
}

static void
check_jx(mfields_base_t *flds)
{
  int bc[3], gdims[3];
  mrc_domain_get_bc(psc.mrc_domain, bc);
  mrc_domain_get_global_dims(psc.mrc_domain, gdims);

  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d_g(p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(p, jx, jy, jz, &ix, &iy, &iz);
      if (ix < 0) {
	if (bc[0] == BC_PERIODIC) {
	  ix += gdims[0];
	} else {
	  continue;
	}
      } else if (ix >= gdims[0]) {
	if (bc[0] == BC_PERIODIC) {
	  ix -= gdims[0];
	} else {
	  continue;
	}
      }
      if (iy < 0) {
	if (bc[1] == BC_PERIODIC) {
	  iy += gdims[1];
	} else {
	  continue;
	}
      } else if (iy >= gdims[1]) {
	if (bc[1] == BC_PERIODIC) {
	  iy -= gdims[1];
	} else {
	  continue;
	}
      }
      if (iz < 0) {
	if (bc[2] == BC_PERIODIC) {
	  iz += gdims[2];
	} else {
	  continue;
	}
      } else if (iz >= gdims[2]) {
	if (bc[2] == BC_PERIODIC) {
	  iz -= gdims[2];
	} else {
	  continue;
	}
      }

      if (F3_BASE(pf, JXI, jx,jy,jz) != iz * 10000 + iy * 100 + ix) {
	printf("ix %d %d %d jx %d %d %d\n", ix, iy, iz, jx, jy, jz);
	printf("exp: %d actual: %g\n", iz * 10000 + iy * 100 + ix,
	       F3_BASE(pf, JXI, jx,jy,jz));
      }
      assert(F3_BASE(pf, JXI, jx,jy,jz) == iz * 10000 + iy * 100 + ix);
    } foreach_3d_end;
  }
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct psc_mod_config conf_c = {
    .mod_bnd = "c",
  };

  // test psc_fill_ghosts()

  psc_create_test_xz(&conf_c);
  mfields_base_t *flds = &psc.flds;
  setup_jx(flds);
  psc_bnd_fill_ghosts(psc.bnd, flds, JXI, JXI + 1);
  check_jx(flds);
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
