
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// mhd_update_finite_volume

static void _mrc_unused
mhd_update_finite_volume(struct ggcm_mhd *mhd, fld3d_t x, fld3d_t fluxes[3],
			 struct mrc_fld *ymask, mrc_fld_data_t dt, int p,
			 int l, int r)
{
  fld3d_foreach(i,j,k, l, r) {
    mrc_fld_data_t ym = ymask ? M3(ymask, 0, i,j,k, p) : 1.f;
    for (int m = 0; m < s_n_comps; m++) {
      F3S(x, m, i,j,k) -= dt * ym *
	(PDE_INV_DX(i) * (F3S(fluxes[0], m, i+di,j,k) - F3S(fluxes[0], m, i,j,k)) +
	 PDE_INV_DY(j) * (F3S(fluxes[1], m, i,j+dj,k) - F3S(fluxes[1], m, i,j,k)) +
	 PDE_INV_DZ(k) * (F3S(fluxes[2], m, i,j,k+dk) - F3S(fluxes[2], m, i,j,k)));
    }
  } fld3d_foreach_end;
}

