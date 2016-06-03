
#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// mhd_update_finite_volume

static void _mrc_unused
mhd_update_finite_volume(struct ggcm_mhd *mhd,
			 struct mrc_fld *x, struct mrc_fld *fluxes[3],
			 struct mrc_fld *ymask, mrc_fld_data_t dt, bool do_correct)
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  if (mhd->amr > 0 && do_correct) {
    for (int d = 0; d < 3; d++) {
      if (gdims[d] > 1) {
	mrc_ddc_amr_apply(mhd->ddc_amr_flux[d], fluxes[d]);
      }
    }
  }

  int nr_comps = mrc_fld_nr_comps(fluxes[0]);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    float *fd1x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, FD1, p);
    float *fd1y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, FD1, p);
    float *fd1z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, FD1, p);

    mrc_fld_foreach(x, i,j,k, 0, 0) {
      mrc_fld_data_t ym = ymask ? M3(ymask, 0, i,j,k, p) : 1.f;
      for (int m = 0; m < nr_comps; m++) {
	M3(x, m, i,j,k, p) -= dt * ym *
	  (fd1x[i] * (M3(fluxes[0], m, i+dx,j,k, p) - M3(fluxes[0], m, i,j,k, p)) +
	   fd1y[j] * (M3(fluxes[1], m, i,j+dy,k, p) - M3(fluxes[1], m, i,j,k, p)) +
	   fd1z[k] * (M3(fluxes[2], m, i,j,k+dz, p) - M3(fluxes[2], m, i,j,k, p)));
      }
    } mrc_fld_foreach_end;
  }
}

