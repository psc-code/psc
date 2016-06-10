
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// mhd_update_finite_volume

static void _mrc_unused
mhd_update_finite_volume(struct ggcm_mhd *mhd,
			 struct mrc_fld *x, struct mrc_fld *fluxes[3],
			 struct mrc_fld *ymask, mrc_fld_data_t dt, bool do_correct,
			 int l, int r)
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  if (do_correct) {
    ggcm_mhd_correct_fluxes(mhd, fluxes);
  }

  int nr_comps = mrc_fld_nr_comps(fluxes[0]);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);

    mrc_fld_foreach(x, i,j,k, l, r) {
      mrc_fld_data_t ym = ymask ? M3(ymask, 0, i,j,k, p) : 1.f;
      for (int m = 0; m < nr_comps; m++) {
	M3(x, m, i,j,k, p) -= dt * ym *
	  (PDE_INV_DX(i) * (M3(fluxes[0], m, i+dx,j,k, p) - M3(fluxes[0], m, i,j,k, p)) +
	   PDE_INV_DY(j) * (M3(fluxes[1], m, i,j+dy,k, p) - M3(fluxes[1], m, i,j,k, p)) +
	   PDE_INV_DZ(k) * (M3(fluxes[2], m, i,j,k+dz, p) - M3(fluxes[2], m, i,j,k, p)));
      }
    } mrc_fld_foreach_end;
  }
}

