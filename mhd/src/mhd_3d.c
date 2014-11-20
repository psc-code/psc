
#include "ggcm_mhd_crds.h"
#include "mhd_util.h"

// ----------------------------------------------------------------------
// compute_B_cc
//
// cell-averaged B

static void __unused
compute_B_cc(struct mrc_fld *B_cc, struct mrc_fld *x, int l, int r)
{
  mrc_fld_foreach(x, i,j,k, l, r) {
    F3(B_cc, 0, i,j,k) = .5f * (F3(x, BX, i,j,k) + F3(x, BX, i+1,j,k));
    F3(B_cc, 1, i,j,k) = .5f * (F3(x, BY, i,j,k) + F3(x, BY, i,j+1,k));
    F3(B_cc, 2, i,j,k) = .5f * (F3(x, BZ, i,j,k) + F3(x, BZ, i,j,k+1));
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// update_ct_uniform

static void __unused
update_ct_uniform(struct ggcm_mhd *mhd,
		  struct mrc_fld *x, struct mrc_fld *E, mrc_fld_data_t dt, int l, int r)
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double dx[3]; mrc_crds_get_dx(crds, dx);
  mrc_fld_data_t dt_on_dx[3] = { dt / dx[0], dt / dx[1], dt / dx[2] };

  const int *ldims = mrc_fld_spatial_dims(x);
  int ie = ldims[0], je = ldims[1], ke = ldims[2];

  for (int k = -l; k < ke + r; k++) {
    for (int j = -l; j < je + r; j++) {
      for (int i = -l; i < ie + r; i++) {
        F3(x, BX, i,j,k) += (dt_on_dx[2] * (F3(E, 1, i  ,j  ,k+1) - F3(E, 1, i  ,j  ,k  )) -
			     dt_on_dx[1] * (F3(E, 2, i  ,j+1,k  ) - F3(E, 2, i  ,j  ,k  )));
        F3(x, BY, i,j,k) += (dt_on_dx[0] * (F3(E, 2, i+1,j  ,k  ) - F3(E, 2, i  ,j  ,k  )) -
			     dt_on_dx[2] * (F3(E, 0, i  ,j  ,k+1) - F3(E, 0, i  ,j  ,k  )));
        F3(x, BZ, i,j,k) += (dt_on_dx[1] * (F3(E, 0, i  ,j+1,k  ) - F3(E, 0, i  ,j  ,k  )) -
			     dt_on_dx[0] * (F3(E, 1, i+1,j  ,k  ) - F3(E, 1, i  ,j  ,k  )));
      }
      int i = ie + r;
      F3(x, BX, i,j,k) += (dt_on_dx[2] * (F3(E, 1, i  ,j  ,k+1) - F3(E, 1, i  ,j  ,k  )) -
			   dt_on_dx[1] * (F3(E, 2, i  ,j+1,k  ) - F3(E, 2, i  ,j  ,k  )));
    }
    for (int i = -l; i < ie + r; i++) {
      int j = je + r;
      F3(x, BY, i,j,k) += (dt_on_dx[0] * (F3(E, 2, i+1,j  ,k  ) - F3(E, 2, i  ,j  ,k  )) -
			   dt_on_dx[2] * (F3(E, 0, i  ,j  ,k+1) - F3(E, 0, i  ,j  ,k  )));
    }
  }
  for (int j = -l; j < je + r; j++) {
    for (int i = -l; i < ie + r; i++) {
      int k = ke + r;
      F3(x, BZ, i,j,k) += (dt_on_dx[1] * (F3(E, 0, i  ,j+1,k  ) - F3(E, 0, i  ,j  ,k  )) -
			   dt_on_dx[0] * (F3(E, 1, i+1,j  ,k  ) - F3(E, 1, i  ,j  ,k  )));
    }
  }
}

// ----------------------------------------------------------------------
// update_ct

static void __unused
update_ct(struct ggcm_mhd *mhd,
	  struct mrc_fld *x, struct mrc_fld *E, mrc_fld_data_t dt)
{
  float *bd3x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bd3y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bd3z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  mrc_fld_foreach(x, i,j,k, 0, 0) {
    F3(x, _B1X, i,j,k) -= dt * (bd3y[j] * (F3(E, 2, i,j+1,k) - F3(E, 2, i,j,k)) -
				bd3z[k] * (F3(E, 1, i,j,k+1) - F3(E, 1, i,j,k)));
    F3(x, _B1Y, i,j,k) -= dt * (bd3z[k] * (F3(E, 0, i,j,k+1) - F3(E, 0, i,j,k)) -
				bd3x[i] * (F3(E, 2, i+1,j,k) - F3(E, 2, i,j,k)));
    F3(x, _B1Z, i,j,k) -= dt * (bd3x[i] * (F3(E, 1, i+1,j,k) - F3(E, 1, i,j,k)) -
				bd3y[j] * (F3(E, 0, i,j+1,k) - F3(E, 0, i,j,k)));
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// mhd_fluxes

static void __unused
mhd_fluxes(struct ggcm_mhd_step *step, struct mrc_fld *fluxes[3], struct mrc_fld *x,
	   struct mrc_fld *B_cc, int bnd, int nghost,
	   void (*flux_func)(struct ggcm_mhd_step *step, struct mrc_fld *fluxes[3],
			     struct mrc_fld *x, struct mrc_fld *B_cc,
			     int ldim, int bnd, int j, int k, int dir))
{
  const int *ldims = mrc_fld_spatial_dims(x);

  for (int k = -bnd; k < ldims[2] + bnd; k++) {
    for (int j = -bnd; j < ldims[1] + bnd; j++) {
      flux_func(step, fluxes, x, B_cc, ldims[0], nghost, j, k, 0);
    }
  }
  for (int k = -bnd; k < ldims[2] + bnd; k++) {
    for (int i = -bnd; i < ldims[0] + bnd; i++) {
      flux_func(step, fluxes, x, B_cc, ldims[1], nghost, k, i, 1);
    }
  }
  for (int j = -bnd; j < ldims[1] + bnd; j++) {
    for (int i = -bnd; i < ldims[0] + bnd; i++) {
      flux_func(step, fluxes, x, B_cc, ldims[2], nghost, i, j, 2);
    }
  }
}

// ----------------------------------------------------------------------
// update_finite_volume_uniform

static void __unused
update_finite_volume_uniform(struct ggcm_mhd *mhd,
			     struct mrc_fld *x, struct mrc_fld *fluxes[3],
			     mrc_fld_data_t dt, int l, int r)
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double dx[3]; mrc_crds_get_dx(crds, dx);
  mrc_fld_data_t dt_on_dx[3] = { dt / dx[0], dt / dx[1], dt / dx[2] };

  mrc_fld_foreach(x, i,j,k, l, r) {
    for (int m = 0; m < 5; m++) {
      F3(x, m, i,j,k) -= 
	(dt_on_dx[0] * (F3(fluxes[0], m, i+1,j,k) - F3(fluxes[0], m, i,j,k)) +
	 dt_on_dx[1] * (F3(fluxes[1], m, i,j+1,k) - F3(fluxes[1], m, i,j,k)) + 
	 dt_on_dx[2] * (F3(fluxes[2], m, i,j,k+1) - F3(fluxes[2], m, i,j,k)));
    }
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// update_finite_volume

static void __unused
update_finite_volume(struct ggcm_mhd *mhd,
		     struct mrc_fld *x, struct mrc_fld *fluxes[3],
		     struct mrc_fld *ymask, mrc_fld_data_t dt)
{
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  mrc_fld_foreach(x, i,j,k, 0, 0) {
    mrc_fld_data_t s = dt * F3(ymask, 0, i,j,k);
    for (int m = 0; m < 5; m++) {
      F3(x, m, i,j,k) -=
	s * (fd1x[i] * (F3(fluxes[0], m, i+1,j,k) - F3(fluxes[0], m, i,j,k)) +
	     fd1y[j] * (F3(fluxes[1], m, i,j+1,k) - F3(fluxes[1], m, i,j,k)) +
	     fd1z[k] * (F3(fluxes[2], m, i,j,k+1) - F3(fluxes[2], m, i,j,k)));
    }
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// mrc_fld_copy_range

// FIXME, mv to right place
static void __unused
mrc_fld_copy_range(struct mrc_fld *to, struct mrc_fld *from, int mb, int me)
{
  assert(to->_nr_ghosts == from->_nr_ghosts);
  int bnd = to->_nr_ghosts;

   mrc_fld_foreach(to, ix,iy,iz, bnd, bnd) {
    for (int m = mb; m < me; m++) {
      F3(to, m, ix,iy,iz) = F3(from, m, ix,iy,iz);
    }
  } mrc_fld_foreach_end;
}

