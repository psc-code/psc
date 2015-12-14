
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_step_private.h"
#include "mhd_util.h"

#include <mrc_ddc.h>

// ----------------------------------------------------------------------
// compute_B_cc
//
// cell-averaged B

static void _mrc_unused
compute_B_cc(struct mrc_fld *B_cc, struct mrc_fld *x, int l, int r)
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    mrc_fld_foreach(x, i,j,k, l, r) {
      M3(B_cc, 0, i,j,k, p) = .5f * (BX_(x, i,j,k, p) + BX_(x, i+dx,j,k, p));
      M3(B_cc, 1, i,j,k, p) = .5f * (BY_(x, i,j,k, p) + BY_(x, i,j+dy,k, p));
      M3(B_cc, 2, i,j,k, p) = .5f * (BZ_(x, i,j,k, p) + BZ_(x, i,j,k+dz, p));
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// compute_Bt_cc
//
// cell-averaged Btotal (ie., add B0 back in, if applicable)

static void _mrc_unused
compute_Bt_cc(struct ggcm_mhd *mhd, struct mrc_fld *B_cc, struct mrc_fld *x, int l, int r)
{
  struct mrc_fld *b0 = mhd->b0;
  if (!b0) {
    compute_B_cc(B_cc, x, l, r);
    return;
  }

  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    mrc_fld_foreach(x, i,j,k, l, r) {
      M3(B_cc, 0, i,j,k, p) = .5f * (BT(x, 0, i,j,k, p) + BT(x, 0, i+dx,j,k, p));
      M3(B_cc, 1, i,j,k, p) = .5f * (BT(x, 1, i,j,k, p) + BT(x, 1, i,j+dy,k, p));
      M3(B_cc, 2, i,j,k, p) = .5f * (BT(x, 2, i,j,k, p) + BT(x, 2, i,j,k+dz, p));
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// update_ct_uniform

void correct_E(struct ggcm_mhd *mhd, struct mrc_fld *E);

static void _mrc_unused
update_ct_uniform(struct ggcm_mhd *mhd,
		  struct mrc_fld *x, struct mrc_fld *E, mrc_fld_data_t dt, int _l, int _r,
		  bool do_correct)
{
  if (mhd->amr > 0 && do_correct) {
    mrc_ddc_amr_apply(mhd->ddc_amr_E, E);
    //    correct_E(mhd, E);
  }

  const float r_db_dt_sq = sqr(mhd->par.r_db_dt);
  int gdims[3]; mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);
  int l[3], r[3];
  for (int d = 0; d < 3; d++) {
    l[d] = (gdims[d] > 1) ? _l : 0;
    r[d] = (gdims[d] > 1) ? _r : 0;
  }

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  const int *ldims = mrc_fld_spatial_dims(x);

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    double ddx[3]; mrc_crds_get_dx(crds, p, ddx);
    mrc_fld_data_t dt_on_dx[3] = { dt / ddx[0], dt / ddx[1], dt / ddx[2] };

    for (int k = -l[2]; k < ldims[2] + r[2]; k++) {
      for (int j = -l[1]; j < ldims[1] + r[1]; j++) {
	for (int i = -l[0]; i < ldims[0] + r[0]; i++) {
	  float crd_fc[3];
	  ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 0, crd_fc);
	  if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	    M3(x, BX, i,j,k, p) += (dt_on_dx[2] * (M3(E, 1, i   ,j   ,k+dz, p) - M3(E, 1, i,j,k, p)) -
				    dt_on_dx[1] * (M3(E, 2, i   ,j+dy,k   , p) - M3(E, 2, i,j,k, p)));
	  }
	  ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 1, crd_fc);
	  if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	    M3(x, BY, i,j,k, p) += (dt_on_dx[0] * (M3(E, 2, i+dx,j   ,k   , p) - M3(E, 2, i,j,k, p)) -
				    dt_on_dx[2] * (M3(E, 0, i   ,j   ,k+dz, p) - M3(E, 0, i,j,k, p)));
	  }
	  ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 2, crd_fc);
	  if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	    M3(x, BZ, i,j,k, p) += (dt_on_dx[1] * (M3(E, 0, i   ,j+dy,k   , p) - M3(E, 0, i,j,k, p)) -
				    dt_on_dx[0] * (M3(E, 1, i+dx,j   ,k   , p) - M3(E, 1, i,j,k, p)));
	  }
	}
        if (gdims[0] > 1) {  
	  int i = ldims[0] + r[0];
	  float crd_fc[3];
	  ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 0, crd_fc);
	  if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	    M3(x, BX, i,j,k, p) += (dt_on_dx[2] * (M3(E, 1, i   ,j   ,k+dz, p) - M3(E, 1, i,j,k, p)) -
				    dt_on_dx[1] * (M3(E, 2, i   ,j+dy,k   , p) - M3(E, 2, i,j,k, p)));
	  }
	}
      }
      if (gdims[1] > 1) {
	for (int i = -l[0]; i < ldims[0] + r[0]; i++) {
	  int j = ldims[1] + r[1];
	  float crd_fc[3];
	  ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 1, crd_fc);
	  if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	    M3(x, BY, i,j,k, p) += (dt_on_dx[0] * (M3(E, 2, i+dx,j   ,k   , p) - M3(E, 2, i,j,k, p)) -
				    dt_on_dx[2] * (M3(E, 0, i   ,j   ,k+dz, p) - M3(E, 0, i,j,k, p)));
	  }
	}
      }
    }
    if (gdims[2] > 1) {
      for (int j = -l[1]; j < ldims[1] + r[1]; j++) {
    	for (int i = -l[0]; i < ldims[0] + r[0]; i++) {
    	  int k = ldims[2] + r[2];
	  float crd_fc[3];
	  ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 2, crd_fc);
	  if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	    M3(x, BZ, i,j,k, p) += (dt_on_dx[1] * (M3(E, 0, i   ,j+dy,k   , p) - M3(E, 0, i,j,k, p)) -
				    dt_on_dx[0] * (M3(E, 1, i+dx,j   ,k   , p) - M3(E, 1, i,j,k, p)));
	  }
    	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// update_ct

static void _mrc_unused
update_ct(struct ggcm_mhd *mhd,
	  struct mrc_fld *x, struct mrc_fld *E, mrc_fld_data_t dt,
	  bool do_correct)
{
  if (mhd->amr > 0 && do_correct) {
    mrc_ddc_amr_apply(mhd->ddc_amr_E, E);
    //    correct_E(mhd, E);
  }

  const float r_db_dt_sq = sqr(mhd->par.r_db_dt);
  int gdims[3]; mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    float *bd3x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, BD3, p);
    float *bd3y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, BD3, p);
    float *bd3z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, BD3, p);

    mrc_fld_foreach(x, i,j,k, 0, 1) {
      float crd_fc[3];
      ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 0, crd_fc);
      if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	M3(x, BX, i,j,k, p) -= dt * (bd3y[j] * (M3(E, 2, i,j+dy,k, p) - M3(E, 2, i,j,k, p)) -
				     bd3z[k] * (M3(E, 1, i,j,k+dz, p) - M3(E, 1, i,j,k, p)));
      }
      ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 1, crd_fc);
      if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	M3(x, BY, i,j,k, p) -= dt * (bd3z[k] * (M3(E, 0, i,j,k+dz, p) - M3(E, 0, i,j,k, p)) -
				     bd3x[i] * (M3(E, 2, i+dx,j,k, p) - M3(E, 2, i,j,k, p)));
      }
      ggcm_mhd_get_crds_fc(mhd, i,j,k, p, 2, crd_fc);
      if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
	M3(x, BZ, i,j,k, p) -= dt * (bd3x[i] * (M3(E, 1, i+dx,j,k, p) - M3(E, 1, i,j,k, p)) -
				     bd3y[j] * (M3(E, 0, i,j+dy,k, p) - M3(E, 0, i,j,k, p)));
      }
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// mhd_fluxes

static void _mrc_unused
mhd_fluxes(struct ggcm_mhd_step *step, struct mrc_fld *fluxes[3], struct mrc_fld *x,
	   struct mrc_fld *B_cc, int bn, int nghost,
	   void (*flux_func)(struct ggcm_mhd_step *step, struct mrc_fld *fluxes[3],
			     struct mrc_fld *x, struct mrc_fld *B_cc,
			     int ldim, int bnd, int j, int k, int dir, int p))
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);

  int bnd[3];
  for (int d = 0; d < 3; d++) {
    if (gdims[d] == 1) {
      bnd[d] = 0;
    } else {
      bnd[d] = bn;
    }
  }

  const int *ldims = mrc_fld_spatial_dims(x);

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    if (gdims[0] > 1) {
      for (int k = -bnd[2]; k < ldims[2] + bnd[2]; k++) {
	for (int j = -bnd[1]; j < ldims[1] + bnd[1]; j++) {
	  flux_func(step, fluxes, x, B_cc, ldims[0], nghost, j, k, 0, p);
	}
      }
    }

    if (gdims[1] > 1) {
      for (int k = -bnd[2]; k < ldims[2] + bnd[2]; k++) {
	for (int i = -bnd[0]; i < ldims[0] + bnd[0]; i++) {
	  flux_func(step, fluxes, x, B_cc, ldims[1], nghost, k, i, 1, p);
	}
      }
    }

    if (gdims[2] > 1) {
      for (int j = -bnd[1]; j < ldims[1] + bnd[1]; j++) {
	for (int i = -bnd[0]; i < ldims[0] + bnd[0]; i++) {
	  flux_func(step, fluxes, x, B_cc, ldims[2], nghost, i, j, 2, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// mhd_fluxes_reconstruct

static void _mrc_unused
mhd_fluxes_reconstruct(struct ggcm_mhd_step *step,
		       struct mrc_fld *U_l[3], struct mrc_fld *U_r[3],
		       struct mrc_fld *x, struct mrc_fld *B_cc, int bn, int nghost,
		       void (*reconstr_func)(struct ggcm_mhd_step *step,
					     struct mrc_fld *U_l[3], struct mrc_fld *U_r[3],
					     struct mrc_fld *x, struct mrc_fld *B_cc,
					     int ldim, int bnd, int j, int k, int dir, int p))
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);

  int bnd[3];
  for (int d = 0; d < 3; d++) {
    if (gdims[d] == 1) {
      bnd[d] = 0;
    } else {
      bnd[d] = bn;
    }
  }

  const int *ldims = mrc_fld_spatial_dims(x);

  // reconstruct
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    if (gdims[0] > 1) {
      for (int k = -bnd[2]; k < ldims[2] + bnd[2]; k++) {
	for (int j = -bnd[1]; j < ldims[1] + bnd[1]; j++) {
	  reconstr_func(step, U_l, U_r, x, B_cc, ldims[0], nghost, j, k, 0, p);
	}
      }
    }

    if (gdims[1] > 1) {
      for (int k = -bnd[2]; k < ldims[2] + bnd[2]; k++) {
	for (int i = -bnd[0]; i < ldims[0] + bnd[0]; i++) {
	  reconstr_func(step, U_l, U_r, x, B_cc, ldims[1], nghost, k, i, 1, p);
	}
      }
    }

    if (gdims[2] > 1) {
      for (int j = -bnd[1]; j < ldims[1] + bnd[1]; j++) {
	for (int i = -bnd[0]; i < ldims[0] + bnd[0]; i++) {
	  reconstr_func(step, U_l, U_r, x, B_cc, ldims[2], nghost, i, j, 2, p);
	}
      }
    }
  }
}

static void _mrc_unused
mhd_fluxes_riemann(struct ggcm_mhd_step *step, struct mrc_fld *fluxes[3],
		   struct mrc_fld *U_l[3], struct mrc_fld *U_r[3],
		   struct mrc_fld *B_cc, int bn, int nghost,
		   void (*riemann_func)(struct ggcm_mhd_step *step, struct mrc_fld *fluxes[3],
					struct mrc_fld *U_l[3], struct mrc_fld *U_r[3],
					struct mrc_fld *B_cc,
					int ldim, int bnd, int j, int k, int dir, int p))
{
  int gdims[3];
  mrc_domain_get_global_dims(fluxes[0]->_domain, gdims);

  int bnd[3];
  for (int d = 0; d < 3; d++) {
    if (gdims[d] == 1) {
      bnd[d] = 0;
    } else {
      bnd[d] = bn;
    }
  }

  const int *ldims = mrc_fld_spatial_dims(fluxes[0]);

  // calculate fluxes
  for (int p = 0; p < mrc_fld_nr_patches(fluxes[0]); p++) {
    if (gdims[0] > 1) {
      for (int k = -bnd[2]; k < ldims[2] + bnd[2]; k++) {
	for (int j = -bnd[1]; j < ldims[1] + bnd[1]; j++) {
	  riemann_func(step, fluxes, U_l, U_r, B_cc, ldims[0], nghost, j, k, 0, p);
	}
      }
    }

    if (gdims[1] > 1) {
      for (int k = -bnd[2]; k < ldims[2] + bnd[2]; k++) {
	for (int i = -bnd[0]; i < ldims[0] + bnd[0]; i++) {
	  riemann_func(step, fluxes, U_l, U_r, B_cc, ldims[1], nghost, k, i, 1, p);
	}
      }
    }

    if (gdims[2] > 1) {
      for (int j = -bnd[1]; j < ldims[1] + bnd[1]; j++) {
	for (int i = -bnd[0]; i < ldims[0] + bnd[0]; i++) {
	  riemann_func(step, fluxes, U_l, U_r, B_cc, ldims[2], nghost, i, j, 2, p);
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// update_finite_volume_uniform

static void _mrc_unused
update_finite_volume_uniform(struct ggcm_mhd *mhd,
			     struct mrc_fld *x, struct mrc_fld *fluxes[3],
			     struct mrc_fld *ymask,
			     mrc_fld_data_t dt, int l, int r, bool do_correct)
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);
  // FIXME ymask!!!

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  if (mhd->amr > 0 && do_correct) {
    for (int d = 0; d < 3; d++) {
      if (gdims[d] > 1) {
	mrc_ddc_amr_apply(mhd->ddc_amr_flux[d], fluxes[d]);
      }
    }
  }

  int nr_comps = mrc_fld_nr_comps(fluxes[0]);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    double ddx[3]; mrc_crds_get_dx(crds, p, ddx);
    mrc_fld_data_t dt_on_dx[3] = { dt / ddx[0], dt / ddx[1], dt / ddx[2] };
    // FIXME, potential for accelerating the 2-d/1-d versions
    for (int d = 0; d < 3; d++) {
      if (gdims[d] == 1) {
	dt_on_dx[d] = 0.f;
      }
    }

    mrc_fld_foreach(x, i,j,k, l, r) {
      mrc_fld_data_t ym = ymask ? M3(ymask, 0, i,j,k, p) : 1.f;
      for (int m = 0; m < nr_comps; m++) {
	M3(x, m, i,j,k, p) -= ym *
	  (dt_on_dx[0] * (M3(fluxes[0], m, i+dx,j,k, p) - M3(fluxes[0], m, i,j,k, p)) +
	   dt_on_dx[1] * (M3(fluxes[1], m, i,j+dy,k, p) - M3(fluxes[1], m, i,j,k, p)) + 
	   dt_on_dx[2] * (M3(fluxes[2], m, i,j,k+dz, p) - M3(fluxes[2], m, i,j,k, p)));
      }
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// update_finite_volume

static void _mrc_unused
update_finite_volume(struct ggcm_mhd *mhd,
		     struct mrc_fld *x, struct mrc_fld *fluxes[3],
		     struct mrc_fld *ymask, mrc_fld_data_t dt, bool do_correct)
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  if (0&&mhd->amr > 0 && do_correct) {
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

// ----------------------------------------------------------------------
// mrc_fld_copy_range

// FIXME, mv to right place
static void _mrc_unused
mrc_fld_copy_range(struct mrc_fld *to, struct mrc_fld *from, int mb, int me)
{
  assert(to->_nr_ghosts == from->_nr_ghosts);
  int bnd = to->_nr_ghosts;

  for (int p = 0; p < mrc_fld_nr_patches(to); p++) {
    mrc_fld_foreach(to, ix,iy,iz, bnd, bnd) {
      for (int m = mb; m < me; m++) {
	M3(to, m, ix,iy,iz, p) = M3(from, m, ix,iy,iz, p);
      }
    } mrc_fld_foreach_end;
  }
}

