
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// mhd_update_finite_volume

static void _mrc_unused
mhd_update_finite_volume(struct ggcm_mhd *mhd,
			 fld3d_t x, struct mrc_fld *fluxes[3],
			 struct mrc_fld *ymask, mrc_fld_data_t dt, int p,
			 int l, int r)
{
  int nr_comps = mrc_fld_nr_comps(fluxes[0]);

  mrc_fld_foreach(ymask, i,j,k, l, r) {
    mrc_fld_data_t ym = ymask ? M3(ymask, 0, i,j,k, p) : 1.f;
    for (int m = 0; m < nr_comps; m++) {
      F3S(x, m, i,j,k) -= dt * ym *
	(PDE_INV_DX(i) * (M3(fluxes[0], m, i+di,j,k, p) - M3(fluxes[0], m, i,j,k, p)) +
	 PDE_INV_DY(j) * (M3(fluxes[1], m, i,j+dj,k, p) - M3(fluxes[1], m, i,j,k, p)) +
	 PDE_INV_DZ(k) * (M3(fluxes[2], m, i,j,k+dk, p) - M3(fluxes[2], m, i,j,k, p)));
    }
  } mrc_fld_foreach_end;
}

