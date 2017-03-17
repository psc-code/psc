
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds.h"

#include <mrc_fld.h>
#include <mrc_domain.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"

#define OPT_FLD1D OPT_FLD1D_C_ARRAY

#include "pde/pde_mhd_divb.c"

static void
ggcm_mhd_calc_divb_cc(struct ggcm_mhd *mhd, struct mrc_fld *f, struct mrc_fld *divB)
{
  fld3d_t p_U, p_divB;
  fld3d_setup(&p_U, f);
  fld3d_setup(&p_divB, divB);

  for (int p = 0; p < mrc_fld_nr_patches(f); p++) {
    pde_patch_set(p);
    fld3d_get(&p_U, p);
    fld3d_get(&p_divB, p);

    patch_calc_divb_bgrid_cc(p_divB, fld3d_make_view(p_U, BX));
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_divb

void
ggcm_mhd_calc_divb(struct ggcm_mhd *mhd, struct mrc_fld *fld, struct mrc_fld *divb)
{
  static bool is_setup = false;
  if (!is_setup) {
    pde_setup(fld);
    pde_mhd_setup(mhd);
    is_setup = true;
  }

  int mhd_type;
  mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);

  struct mrc_fld *f = mrc_fld_get_as(fld, FLD_TYPE);
  struct mrc_fld *d = mrc_fld_get_as(divb, FLD_TYPE);

  if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
    ggcm_mhd_calc_divb_cc(mhd, f, d);
    goto out;
  }
  
  mrc_fld_data_t hx = 1., hy = 1., hz = 1.;
  int gdims[3];
  mrc_domain_get_global_dims(fld->_domain, gdims);
  if (gdims[0] == 1) hx = 0.;
  if (gdims[1] == 1) hy = 0.;
  if (gdims[2] == 1) hz = 0.;
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);


  struct mrc_fld *ymask = NULL;
  if (mhd->ymask) {
    ymask = mrc_fld_get_as(mhd->ymask, FLD_TYPE);
  }

  mrc_fld_data_t max = 0.;

  if (MT_BGRID(mhd_type) == MT_BGRID_FC_GGCM) {
    for (int p = 0; p < mrc_fld_nr_patches(divb); p++) {
      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(fld->_domain, p, &info);
  
      float *bd3x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD3, p);
      float *bd3y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD3, p);
      float *bd3z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD3, p);

      mrc_fld_foreach(divb, ix,iy,iz, 0, 0) {
	M3(d,0, ix,iy,iz, p) =
	  (BX_(f, ix,iy,iz, p) - BX_(f, ix-dx,iy,iz, p)) * bd3x[ix] +
	  (BY_(f, ix,iy,iz, p) - BY_(f, ix,iy-dy,iz, p)) * bd3y[iy] +
	  (BZ_(f, ix,iy,iz, p) - BZ_(f, ix,iy,iz-dz, p)) * bd3z[iz];
 	
	// the incoming solar wind won't match and hence divb != 0 here
	if ((info.off[0] == 0 && ix <= 0) || M3(ymask, 0, ix,iy,iz, p) == 0)
	  continue;
	max = mrc_fld_max(max, mrc_fld_abs(M3(d,0, ix,iy,iz, p)));
      } mrc_fld_foreach_end;
    }
  } else if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
    for (int p = 0; p < mrc_fld_nr_patches(divb); p++) {
      float *bd3x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD3, p);
      float *bd3y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD3, p);
      float *bd3z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD3, p);

      mrc_fld_foreach(divb, ix,iy,iz, 0, 0) {
	M3(d,0, ix,iy,iz, p) =
	  (BX_(f, ix+dx,iy,iz, p) - BX_(f, ix,iy,iz, p)) * hx * bd3x[ix] +
	  (BY_(f, ix,iy+dy,iz, p) - BY_(f, ix,iy,iz, p)) * hy * bd3y[iy] +
	  (BZ_(f, ix,iy,iz+dz, p) - BZ_(f, ix,iy,iz, p)) * hz * bd3z[iz];

	if (M3(ymask, 0, ix,iy,iz, p) == 0)
	  continue;
	max = mrc_fld_max(max, mrc_fld_abs(M3(d,0, ix,iy,iz, p)));
      } mrc_fld_foreach_end;
    }
  } else {
    assert(0);
  }

  if (mhd->ymask) {
    mrc_fld_put_as(ymask, mhd->ymask);
  }

  mrc_fld_data_t max_divb;
  MPI_Allreduce(&max, &max_divb, 1, MPI_MRC_FLD_DATA_T, MPI_MAX, ggcm_mhd_comm(mhd));
  mpi_printf(ggcm_mhd_comm(mhd), "max divb = %g\n", max_divb);

 out:
  mrc_fld_put_as(f, fld);
  mrc_fld_put_as(d, divb);

  if (0) {
    pde_free();
  }
}

