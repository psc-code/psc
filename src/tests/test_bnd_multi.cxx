
#include "psc_testing.h"
#include "psc_bnd.h"
#include "psc_fields_as_c.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#if 0
static void
setup_jx(struct psc_mfields *mflds_base)
{
  mfields_t mf = mflds_base->get_as<mfields_t>(0, 0);
  psc_foreach_patch(ppsc, p) {
    Fields F(mf[p]);
    psc_foreach_3d(ppsc, p, jx, jy, jz, 0, 0) {
      int ix, iy, iz;
      psc_local_to_global_indices(ppsc, p, jx, jy, jz, &ix, &iy, &iz);
      F(JXI, jx,jy,jz) = iz * 10000 + iy * 100 + ix;
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, JXI, JXI + 1);
}

static void
check_jx(struct psc_mfields *mflds_base)
{
  int bc[3], gdims[3];
  mrc_domain_get_bc(ppsc->mrc_domain, bc);
  mrc_domain_get_global_dims(ppsc->mrc_domain, gdims);

  mfields_t mf = mflds_base->get_as<mfields_t>(0, 0);
  psc_foreach_patch(ppsc, p) {
    Fields F(mf[p]);
    psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
      int ix, iy, iz;
      psc_local_to_global_indices(ppsc, p, jx, jy, jz, &ix, &iy, &iz);
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

      if (F(JXI, jx,jy,jz) != iz * 10000 + iy * 100 + ix) {
	printf("ix %d %d %d jx %d %d %d\n", ix, iy, iz, jx, jy, jz);
	printf("exp: %d actual: %g\n", iz * 10000 + iy * 100 + ix,
	       F(JXI, jx,jy,jz));
      }
      assert(F(JXI, jx,jy,jz) == iz * 10000 + iy * 100 + ix);
    } foreach_3d_end;
  }
  mf.put_as(mflds_base, 0, 0);
}
#endif

int
main(int argc, char **argv)
{
#if 0
  psc_testing_init(&argc, &argv);

  // test psc_fill_ghosts()

  struct psc_case *_case = psc_create_test_xz();
  psc_bnd_set_type(ppsc->bnd, "c");
  psc_case_setup(_case);
  struct psc_mfields *flds = ppsc->flds;
  setup_jx(flds);
  psc_bnd_fill_ghosts(ppsc->bnd, flds, JXI, JXI + 1);
  check_jx(flds);
  psc_case_destroy(_case);

  psc_testing_finalize();
#endif
}
