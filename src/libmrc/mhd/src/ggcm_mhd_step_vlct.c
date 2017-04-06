
#include <ggcm_mhd_step_private.h>

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <ggcm_mhd_diag.h>

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_io.h>

#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"

#define OPT_EQN OPT_EQN_MHD_FCONS
#define OPT_BACKGROUND false
#define OPT_GET_DT OPT_GET_DT_MHD_CT

// FIXME, temp hack that should be done as not compile-time regular
// option some time
#define OPT_DIVB_CT

#include "pde/pde_setup.c"
#include "pde/pde_mhd_setup.c"
#include "pde/pde_mhd_line.c"
#include "pde/pde_mhd_convert.c"
#include "pde/pde_mhd_reconstruct.c"
#include "pde/pde_mhd_divb_glm.c"
#include "pde/pde_mhd_riemann.c"
#include "pde/pde_mhd_stage.c"
#include "pde/pde_mhd_get_dt.c"

#include "mhd_3d.c"


static const int nghost = 4; // FIXME

#include <string.h>

static int ldims[3];

// ======================================================================
// ggcm_mhd_step subclass "vlct"

struct ggcm_mhd_step_vlct {
  struct mhd_options opt;

  bool debug_dump;

  fld1d_state_t U;
  fld1d_state_t U_l;
  fld1d_state_t U_r;
  fld1d_state_t W;
  fld1d_state_t W_l;
  fld1d_state_t W_r;
  fld1d_state_t F;
  fld1d_t bx;
};

#define ggcm_mhd_step_vlct(step) mrc_to_subobj(step, struct ggcm_mhd_step_vlct)

// ======================================================================

static inline mrc_fld_data_t
l_r_avg(mrc_fld_data_t u, mrc_fld_data_t l, mrc_fld_data_t r)
{
  if (u > 0.) {
    return l;
  } else if (u < 0.) {
    return r;
  } else {
    return .5 * (l + r);
  }
}

#define MAKE_COMPUTE_E_EDGE(X,Y,Z, _IJ,_JJ,_KJ, _IK,_JK,_KK)		\
									\
static void								\
compute_E##X##_edge(struct mrc_fld *E, struct mrc_fld *Ecc,		\
		    struct mrc_fld **fluxes, int bnd)			\
{									\
  int gdims[3];								\
  mrc_domain_get_global_dims(E->_domain, gdims);			\
									\
  int BY = BX + Y, BZ = BX + Z;						\
  int IJ = gdims[0] > 1 ? _IJ : 0;					\
  int JJ = gdims[1] > 1 ? _JJ : 0;					\
  int KJ = gdims[2] > 1 ? _KJ : 0;					\
  int IK = gdims[0] > 1 ? _IK : 0;					\
  int JK = gdims[1] > 1 ? _JK : 0;					\
  int KK = gdims[2] > 1 ? _KK : 0;					\
									\
  for (int p = 0; p < mrc_fld_nr_patches(E); p++) {			\
    mrc_fld_foreach(E, i,j,k, bnd, bnd+1) {				\
      mrc_fld_data_t de1_l2, de1_r2, de1_l3, de1_r3;			\
      									\
      de1_l3 = l_r_avg(M3(fluxes[Y], RR, i-IK,j-JK,k-KK, p),		\
		       M3(fluxes[Z], BY, i   ,j   ,k   , p) - M3(Ecc, X, i   -IK,j   -JK,k   -KK, p), \
		       M3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ, p) - M3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK, p)); \
      									\
      de1_r3 = l_r_avg(M3(fluxes[Y], RR, i,j,k, p),			\
		       M3(fluxes[Z], BY, i   ,j   ,k   , p) - M3(Ecc, X, i   ,j   ,k   , p), \
		       M3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ, p) - M3(Ecc, X, i-IJ,j-JJ,k-KJ, p)); \
      									\
      de1_l2 = l_r_avg( M3(fluxes[Z], RR, i-IJ,j-JJ,k-KJ, p),		\
		       -M3(fluxes[Y], BZ, i   ,j   ,k   , p) - M3(Ecc, X, i-IJ   ,j-JJ   ,k-KJ   , p), \
		       -M3(fluxes[Y], BZ, i-IK,j-JK,k-KK, p) - M3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK, p)); \
      									\
      de1_r2 = l_r_avg( M3(fluxes[Z], RR, i,j,k, p),			\
		       -M3(fluxes[Y], BZ, i   ,j   ,k   , p) - M3(Ecc, X, i   ,j   ,k   , p), \
		       -M3(fluxes[Y], BZ, i-IK,j-JK,k-KK, p) - M3(Ecc, X, i-IK,j-JK,k-KK, p)); \
      									\
      M3(E, X, i,j,k, p) =						\
	.25 * (  M3(fluxes[Z], BY, i,j,k, p) + M3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ, p) \
	       - M3(fluxes[Y], BZ, i,j,k, p) - M3(fluxes[Y], BZ, i-IK,j-JK,k-KK, p) \
	       + de1_l2 + de1_r2 + de1_l3 + de1_r3);		        \
    } mrc_fld_foreach_end;						\
  }									\
}

#define MAKE_COMPUTE_E_EDGE_(X,Y,Z, _IJ,_JJ,_KJ, _IK,_JK,_KK)		\
									\
static void								\
compute_E##X##_edge(struct mrc_fld *E, struct mrc_fld *Ecc,		\
		   struct mrc_fld **fluxes, int bnd)			\
{									\
  int gdims[3];								\
  mrc_domain_get_global_dims(E->_domain, gdims);			\
									\
  int BY = BX + Y, BZ = BX + Z;						\
  int IJ = gdims[0] > 1 ? _IJ : 0;					\
  int JJ = gdims[1] > 1 ? _JJ : 0;					\
  int KJ = gdims[2] > 1 ? _KJ : 0;					\
  int IK = gdims[0] > 1 ? _IK : 0;					\
  int JK = gdims[1] > 1 ? _JK : 0;					\
  int KK = gdims[2] > 1 ? _KK : 0;					\
									\
  for (int p = 0; p < mrc_fld_nr_patches(E); p++) {			\
    mrc_fld_foreach(E, i,j,k, bnd, bnd+1) {				\
      mrc_fld_data_t de1_l2, de1_r2, de1_l3, de1_r3;			\
      									\
      de1_l3 = l_r_avg(M3(fluxes[Y], RR, i-IK,j-JK,k-KK, p),		\
		       M3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ, p) - M3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK, p), \
		       M3(fluxes[Z], BY, i   ,j   ,k   , p) - M3(Ecc, X, i   -IK,j   -JK,k   -KK, p)); \
      									\
      de1_r3 = l_r_avg(M3(fluxes[Y], RR, i,j,k, p),			\
		       M3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ, p) - M3(Ecc, X, i-IJ,j-JJ,k-KJ, p), \
		       M3(fluxes[Z], BY, i   ,j   ,k   , p) - M3(Ecc, X, i   ,j   ,k   , p)); \
      									\
      de1_l2 = l_r_avg( M3(fluxes[Z], RR, i-IJ,j-JJ,k-KJ, p),		\
		       -M3(fluxes[Y], BZ, i-IK,j-JK,k-KK, p) - M3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK, p), \
		       -M3(fluxes[Y], BZ, i   ,j   ,k   , p) - M3(Ecc, X, i-IJ   ,j-JJ   ,k-KJ   , p)); \
      									\
      de1_r2 = l_r_avg( M3(fluxes[Z], RR, i,j,k, p),			\
		       -M3(fluxes[Y], BZ, i-IK,j-JK,k-KK, p) - M3(Ecc, X, i-IK,j-JK,k-KK, p), \
		       -M3(fluxes[Y], BZ, i   ,j   ,k   , p) - M3(Ecc, X, i   ,j   ,k   , p)); \
    									\
      M3(E, X, i,j,k, p) =						\
	.25 * (  M3(fluxes[Z], BY, i,j,k, p) + M3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ, p) \
	       - M3(fluxes[Y], BZ, i,j,k, p) - M3(fluxes[Y], BZ, i-IK,j-JK,k-KK, p) \
	       + de1_l2 + de1_r2 + de1_l3 + de1_r3);		        \
    } mrc_fld_foreach_end;						\
  }									\
}

MAKE_COMPUTE_E_EDGE_(0,1,2, 0,1,0, 0,0,1)
MAKE_COMPUTE_E_EDGE_(1,2,0, 0,0,1, 1,0,0)
MAKE_COMPUTE_E_EDGE_(2,0,1, 1,0,0, 0,1,0)

static void
compute_E(struct ggcm_mhd_step *step, struct mrc_fld *E,
	  struct mrc_fld *x, struct mrc_fld *B_cc, struct mrc_fld *fluxes[3],
	  int bnd)
{
  struct mrc_fld *Ecc = ggcm_mhd_get_3d_fld(step->mhd, 3);

  // calculate cell-centered E first
  for (int p = 0; p < mrc_fld_nr_patches(Ecc); p++) {
    mrc_fld_foreach(Ecc, i,j,k, bnd, bnd) {
      mrc_fld_data_t rri = 1. / RR_(x, i,j,k, p);
      mrc_fld_data_t B[3] = { M3(B_cc, 0, i,j,k, p), M3(B_cc, 1, i,j,k, p), M3(B_cc, 2, i,j,k, p) };
      mrc_fld_data_t v[3] = { rri * RVX_(x, i,j,k, p), rri * RVY_(x, i,j,k, p), rri * RVZ_(x, i,j,k, p) };
      M3(Ecc, 0, i,j,k, p) = B[1] * v[2] - B[2] * v[1];
      M3(Ecc, 1, i,j,k, p) = B[2] * v[0] - B[0] * v[2];
      M3(Ecc, 2, i,j,k, p) = B[0] * v[1] - B[1] * v[0];
    } mrc_fld_foreach_end;
  }

  // then calculate edge-centered E based on the cell-centered one and
  // the fluxes
  compute_E0_edge(E, Ecc, fluxes, bnd - 1);
  compute_E1_edge(E, Ecc, fluxes, bnd - 1);
  compute_E2_edge(E, Ecc, fluxes, bnd - 1);

  ggcm_mhd_put_3d_fld(step->mhd, Ecc);
}

static void
fluxes_pred(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  fld1d_state_t U = sub->U, U_l = sub->U_l, U_r = sub->U_r;
  fld1d_state_t W = sub->W, W_l = sub->W_l, W_r = sub->W_r;
  fld1d_state_t F = sub->F;
  fld1d_t bx = sub->bx;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_for_each_dir(dir) {
      int ldim = s_ldims[dir];
      pde_for_each_line(dir, j, k, nghost) {
	mhd_get_line_state_fcons_ct(U, bx, x, B_cc, j, k, dir, p, -nghost, ldim + nghost);
	mhd_prim_from_cons(W, U, -nghost, ldim + nghost);
	mhd_reconstruct_pcm(U_l, U_r, W_l, W_r, W, bx, -(nghost - 1), ldim + nghost);
	mhd_riemann(F, U_l, U_r, W_l, W_r, -(nghost - 1), ldim + nghost);
	mhd_put_line_state_fcons_ct(flux[dir], F, j, k, dir, p, -(nghost - 1), ldim + nghost);
      }
    }
  }
}

static void
fluxes_corr(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  fld1d_state_t U = sub->U, U_l = sub->U_l, U_r = sub->U_r;
  fld1d_state_t W = sub->W, W_l = sub->W_l, W_r = sub->W_r;
  fld1d_state_t F = sub->F;
  fld1d_t bx = sub->bx;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_for_each_dir(dir) {
      int ldim = s_ldims[dir];
      pde_for_each_line(dir, j, k, 1) {
	mhd_get_line_state_fcons_ct(U, bx, x, B_cc, j, k, dir, p, - (nghost - 1), ldim + (nghost - 1));
	mhd_prim_from_cons(W, U, - (nghost - 1), ldim + (nghost - 1));
	mhd_reconstruct(U_l, U_r, W_l, W_r, W, bx, 0, ldim + 1);
	mhd_riemann(F, U_l, U_r, W_l, W_r, 0, ldim + 1);
	mhd_put_line_state_fcons_ct(flux[dir], F, j, k, dir, p, 0, ldim + 1);
      }
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_setup

static void
ggcm_mhd_step_vlct_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));

  fld1d_state_setup(&sub->U);
  fld1d_state_setup(&sub->U_l);
  fld1d_state_setup(&sub->U_r);
  fld1d_state_setup(&sub->W);
  fld1d_state_setup(&sub->W_l);
  fld1d_state_setup(&sub->W_r);
  fld1d_state_setup(&sub->F);
  fld1d_setup(&sub->bx);

  mhd->ymask = ggcm_mhd_get_3d_fld(mhd, 1);
  mrc_fld_set(mhd->ymask, 1.);

  ggcm_mhd_step_setup_member_objs_sub(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_destroy

static void
ggcm_mhd_step_vlct_destroy(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_put_3d_fld(mhd, mhd->ymask);

  pde_free();
}

// ----------------------------------------------------------------------
// compute_Ediffu_const

static void
compute_Ediffu_const(struct ggcm_mhd_step *step, struct mrc_fld *E_ec,
                     struct mrc_fld *x, struct mrc_fld *B_cc)
{
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *j_ec = ggcm_mhd_get_3d_fld(mhd, 3);
  struct mrc_fld *j_cc = ggcm_mhd_get_3d_fld(mhd, 3);
  mrc_fld_data_t Jx_ecy, Jx_ecz, Jy_ecx, Jy_ecz, Jz_ecx, Jz_ecy;
  mrc_fld_data_t Bx_ecy, Bx_ecz, By_ecx, By_ecz, Bz_ecx, Bz_ecy;

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int II = gdims[0] > 1 ? 1 : 0;
  int JJ = gdims[1] > 1 ? 1 : 0;
  int KK = gdims[2] > 1 ? 1 : 0;

  // calc edge-centered J
  for (int p = 0; p < mrc_fld_nr_patches(j_ec); p++) {
    double dx[3]; mrc_crds_get_dx(crds, p, dx);
    mrc_fld_data_t dxi[3] = { 1. / dx[0], 1. / dx[1], 1. / dx[2] };
    // don't divide by 0 in invariant dimensions

    mrc_fld_foreach(j_ec, i,j,k, 1, 2) {
      M3(j_ec, 0, i,j,k, p) =
	(BZ_(x, i,j,k, p) - BZ_(x, i   , j-JJ, k   , p)) * dxi[1] -
	(BY_(x, i,j,k, p) - BY_(x, i   , j   , k-KK, p)) * dxi[2];
      M3(j_ec, 1, i,j,k, p) =
	(BX_(x, i,j,k, p) - BX_(x, i   , j   , k-KK, p)) * dxi[2] -
	(BZ_(x, i,j,k, p) - BZ_(x, i-II, j   , k   , p)) * dxi[0];
      M3(j_ec, 2, i,j,k, p) =
	(BY_(x, i,j,k, p) - BY_(x, i-II, j   , k   , p)) * dxi[0] -
	(BX_(x, i,j,k, p) - BX_(x, i   , j-JJ, k   , p)) * dxi[1];
    } mrc_fld_foreach_end;
  }

  // calc cell-centered J
  for (int p = 0; p < mrc_fld_nr_patches(E_ec); p++) {
    mrc_fld_foreach(j_ec, i,j,k, 1, 1) {
      M3(j_cc, 0, i,j,k, p) =
	0.25 * (M3(j_ec, 0, i   , j+JJ, k+KK, p) + M3(j_ec, 0, i   , j   , k+KK, p) +
		M3(j_ec, 0, i   , j+JJ, k   , p) + M3(j_ec, 0, i   , j   , k   , p));
      M3(j_cc, 1, i,j,k, p) =
	0.25 * (M3(j_ec, 1, i+II, j   , k+KK, p) + M3(j_ec, 1, i+II, j   , k   , p) +
		M3(j_ec, 1, i   , j   , k+KK, p) + M3(j_ec, 1, i   , j   , k   , p));
      M3(j_cc, 2, i,j,k, p) =
	0.25 * (M3(j_ec, 2, i+II, j+JJ, k   , p) + M3(j_ec, 2, i+II, j   , k   , p) +
		M3(j_ec, 2, i   , j+JJ, k   , p) + M3(j_ec, 2, i   , j   , k   , p));
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t eta = mhd->par.diffco / mhd->resnorm;
  mrc_fld_data_t d_i = mhd->par.d_i;

  for (int p = 0; p < mrc_fld_nr_patches(E_ec); p++) {
    mrc_fld_foreach(j_ec, i,j,k, 0, 1) {
      // eta J
      M3(E_ec, 0, i,j,k, p) = eta * M3(j_ec, 0, i,j,k, p);
      M3(E_ec, 1, i,j,k, p) = eta * M3(j_ec, 1, i,j,k, p);
      M3(E_ec, 2, i,j,k, p) = eta * M3(j_ec, 2, i,j,k, p);

      // d_i J x B
      if (d_i > 0.0) {
	// average edge centered J to the edges needed for JxB
	// the ec_[xyz] says which edge J is on, aka the component of E that
	// the value is used to calculate
	// Jy_ecx = 0.25 * (M3(j_ec, 1, i+II, j-JJ, k   , p) + M3(j_ec, 1, i+II, j   , k  , p) +
	// 		    M3(j_ec, 1, i   , j-JJ, k   , p) + M3(j_ec, 1, i   , j   , k  , p));
	// Jz_ecx = 0.25 * (M3(j_ec, 2, i+II, j   , k-KK, p) + M3(j_ec, 2, i+II, j   , k  , p) +
	// 		    M3(j_ec, 2, i   , j   , k-KK, p) + M3(j_ec, 2, i   , j   , k  , p));
	// Jx_ecy = 0.25 * (M3(j_ec, 0, i-II, j+JJ, k   , p) + M3(j_ec, 0, i-II, j   , k  , p) +
	// 		    M3(j_ec, 0, i   , j+JJ, k   , p) + M3(j_ec, 0, i   , j   , k  , p));
	// Jz_ecy = 0.25 * (M3(j_ec, 2, i   , j+JJ, k-KK, p) + M3(j_ec, 2, i   , j+JJ, k  , p) +
	// 		    M3(j_ec, 2, i   , j   , k-KK, p) + M3(j_ec, 2, i   , j   , k  , p));
	// Jx_ecz = 0.25 * (M3(j_ec, 0, i-II, j   , k+KK, p) + M3(j_ec, 0, i-II, j   , k  , p) +
	// 		    M3(j_ec, 0, i   , j   , k+KK, p) + M3(j_ec, 0, i   , j   , k  , p));
	// Jy_ecz = 0.25 * (M3(j_ec, 1, i   , j-JJ, k+KK, p) + M3(j_ec, 1, i   , j-JJ, k  , p) +
	// 		    M3(j_ec, 1, i   , j   , k+KK, p) + M3(j_ec, 1, i   , j   , k  , p));

	// use j_cc to get J used for j x b
	Jy_ecx = 0.25 * (M3(j_cc, 1, i   , j-JJ, k-KK, p) + M3(j_cc, 1, i   , j-JJ, k   , p) +
			 M3(j_cc, 1, i   , j   , k-KK, p) + M3(j_cc, 1, i   , j   , k   , p));
	Jz_ecx = 0.25 * (M3(j_cc, 2, i   , j-JJ, k-KK, p) + M3(j_cc, 2, i   , j-JJ, k   , p) +
			 M3(j_cc, 2, i   , j   , k-KK, p) + M3(j_cc, 2, i   , j   , k   , p));
	Jx_ecy = 0.25 * (M3(j_cc, 0, i-II, j   , k-KK, p) + M3(j_cc, 0, i-II, j   , k   , p) +
			 M3(j_cc, 0, i-II, j   , k   , p) + M3(j_cc, 0, i   , j   , k   , p));
	Jz_ecy = 0.25 * (M3(j_cc, 2, i-II, j   , k-KK, p) + M3(j_cc, 2, i-II, j   , k   , p) +
			 M3(j_cc, 2, i-II, j   , k   , p) + M3(j_cc, 2, i   , j   , k   , p));
	Jx_ecz = 0.25 * (M3(j_cc, 0, i-II, j-JJ, k   , p) + M3(j_cc, 0, i-II, j   , k   , p) +
			 M3(j_cc, 0, i-II, j   , k   , p) + M3(j_cc, 0, i   , j   , k   , p));
	Jy_ecz = 0.25 * (M3(j_cc, 1, i-II, j-JJ, k   , p) + M3(j_cc, 1, i-II, j   , k   , p) +
			 M3(j_cc, 1, i-II, j   , k   , p) + M3(j_cc, 1, i   , j   , k   , p));

	// average face centered B to edge centers
	By_ecx = 0.5 * (BY_(x, i   , j   , k-KK, p) + BY_(x, i, j, k, p));
	Bz_ecx = 0.5 * (BZ_(x, i   , j-JJ, k   , p) + BZ_(x, i, j, k, p));
	Bx_ecy = 0.5 * (BX_(x, i   , j   , k-KK, p) + BX_(x, i, j, k, p));
	Bz_ecy = 0.5 * (BZ_(x, i-II, j   , k   , p) + BZ_(x, i, j, k, p));
	Bx_ecz = 0.5 * (BX_(x, i   , j-JJ, k   , p) + BX_(x, i, j, k, p));
	By_ecz = 0.5 * (BY_(x, i-II, j   , k   , p) + BY_(x, i, j, k, p));
	
	M3(E_ec, 0, i,j,k, p) += d_i * ( Jy_ecx * Bz_ecx - Jz_ecx * By_ecx);
	M3(E_ec, 1, i,j,k, p) += d_i * (-Jx_ecy * Bz_ecy + Jz_ecy * Bx_ecy);
	M3(E_ec, 2, i,j,k, p) += d_i * ( Jx_ecz * By_ecz - Jy_ecz * Bx_ecz);
      }
    } mrc_fld_foreach_end;
  }

  // FIXME!!! energy flux

  ggcm_mhd_put_3d_fld(mhd, j_ec);
  ggcm_mhd_put_3d_fld(mhd, j_cc);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_get_dt

static double
ggcm_mhd_step_vlct_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  return pde_mhd_get_dt(step->mhd, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_run

static void
ggcm_mhd_step_vlct_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *B_cc = ggcm_mhd_get_3d_fld(mhd, 3);

  ldims[0] = mrc_fld_spatial_dims(x)[0];
  ldims[1] = mrc_fld_spatial_dims(x)[1];
  ldims[2] = mrc_fld_spatial_dims(x)[2];

  struct mrc_fld *x_half = ggcm_mhd_get_3d_fld(mhd, 8);
  struct mrc_fld *E_ec = ggcm_mhd_get_3d_fld(mhd, 3);
  struct mrc_fld *flux[3] = { ggcm_mhd_get_3d_fld(mhd, 8),
			      ggcm_mhd_get_3d_fld(mhd, 8),
			      ggcm_mhd_get_3d_fld(mhd, 8), };

  // FIXME, this is done in get_dt, and redoing it could be avoided
  ggcm_mhd_fill_ghosts(mhd, x, mhd->time_code);

  mrc_fld_data_t dt = mhd->dt_code;

  // resistivity

  // FIXME, this is done in get_dt, and redoing it could be avoided
  compute_B_cc(B_cc, x, 3, 3);
  if (mhd->par.magdiffu == MAGDIFFU_CONST) {
    if (mhd->par.diffco > 0.) {
      compute_Ediffu_const(step, E_ec, x, B_cc);
      ggcm_mhd_fill_ghosts_E(mhd, E_ec);
      update_ct_uniform(mhd, x, E_ec, dt, 0, 0, false);
      ggcm_mhd_fill_ghosts(mhd, x, mhd->time_code);
      compute_B_cc(B_cc, x, 3, 3);
    }
  } else {
    mpi_printf(ggcm_mhd_comm(mhd), "WARNING: magdiffu '%d' not handled!\n",
	       mhd->par.magdiffu);
  }

  // PREDICTOR

  if (sub->debug_dump) {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_type(diag, "c");
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "dbg");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1:rr:v:pp:b:divb");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_fill_ghosts(mhd, x, mhd->time_code);
    ggcm_mhd_diag_run_now(diag, x, DIAG_TYPE_3D, cnt++);
  }

  fld3d_t _x_half, _x, ymask, _flux[3];
  fld3d_setup(&_x_half, x_half);
  fld3d_setup(&_x, x);
  fld3d_setup(&ymask, mhd->ymask);
  for (int d = 0; d < 3; d++) {
    fld3d_setup(&_flux[d], flux[d]);
  }

  mrc_fld_copy_range(x_half, x, 0, 8);
  fluxes_pred(step, flux, x, B_cc);
  for (int p = 0; p < mrc_fld_nr_patches(x_half); p++) {
    pde_patch_set(p);
    fld3d_get(&_x_half, p);
    fld3d_get(&ymask, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&_flux[d], p);
    }

    mhd_update_finite_volume(mhd, _x_half, _flux, ymask, .5 * dt, nghost - 1, nghost - 1);

    fld3d_put(&_x_half, p);
    fld3d_put(&ymask, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&_flux[d], p);
    }
  }
  compute_E(step, E_ec, x, B_cc, flux, 4);
  ggcm_mhd_fill_ghosts_E(mhd, E_ec);
  update_ct_uniform(mhd, x_half, E_ec, .5 * dt, nghost - 1, nghost - 1, false);

  // CORRECTOR

  compute_B_cc(B_cc, x_half, nghost - 1, nghost - 1);
  fluxes_corr(step, flux, x_half, B_cc);
  ggcm_mhd_correct_fluxes(mhd, flux);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);
    fld3d_get(&_x, p);
    fld3d_get(&ymask, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&_flux[d], p);
    }

    mhd_update_finite_volume(mhd, _x, _flux, ymask, dt, 0, 0);

    fld3d_put(&_x, p);
    fld3d_put(&ymask, p);
    for (int d = 0; d < 3; d++) {
      fld3d_put(&_flux[d], p);
    }
  }
  compute_E(step, E_ec, x_half, B_cc, flux, 4);
  ggcm_mhd_fill_ghosts_E(mhd, E_ec);
  update_ct_uniform(mhd, x, E_ec, dt, 0, 0, true);

  // clean up

  ggcm_mhd_put_3d_fld(mhd, B_cc);
  ggcm_mhd_put_3d_fld(mhd, x_half);
  ggcm_mhd_put_3d_fld(mhd, E_ec);
  ggcm_mhd_put_3d_fld(mhd, flux[0]);
  ggcm_mhd_put_3d_fld(mhd, flux[1]);
  ggcm_mhd_put_3d_fld(mhd, flux[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_setup_flds

static void
ggcm_mhd_step_vlct_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 4);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_FCONS_FC);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_get_e_ec
//
// This is very heavy for for just calculating E, and it's not strictly
// speaking the same E as used in the time step due to the operator
// splitting.

static void
ggcm_mhd_step_vlct_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *Eout,
                            struct mrc_fld *state_vec)
{
  // the state vector should already be FLD_TYPE, but Eout is the data type
  // of the output
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *E = mrc_fld_get_as(Eout, FLD_TYPE);
  struct mrc_fld *x = mrc_fld_get_as(state_vec, FLD_TYPE);

  struct mrc_fld *B_cc = ggcm_mhd_get_3d_fld(mhd, 3);
  struct mrc_fld *Ediff = ggcm_mhd_get_3d_fld(mhd, 3);
  struct mrc_fld *Econv = ggcm_mhd_get_3d_fld(mhd, 3);
  struct mrc_fld *flux[3] = { ggcm_mhd_get_3d_fld(mhd, 8),
                              ggcm_mhd_get_3d_fld(mhd, 8),
                              ggcm_mhd_get_3d_fld(mhd, 8), };

  // Do diffusive terms first
  compute_B_cc(B_cc, x, 3, 3);
  if (mhd->par.magdiffu == MAGDIFFU_CONST) {
    compute_Ediffu_const(step, Ediff, x, B_cc);
  } else {
    for (int p = 0; p < mrc_fld_nr_patches(Ediff); p++) {
      mrc_fld_foreach(Ediff, i,j,k, 2, 2) {
	for (int d=0; d < 3; d++) {
	  M3(Ediff, d, i,j,k, p) = 0.0;
	}
      }  mrc_fld_foreach_end;
    }
    mpi_printf(ggcm_mhd_comm(mhd), "WARNING: magdiffu '%d' not handled for diags!\n",
	       mhd->par.magdiffu);
  }

  // Do convective term using the apprepriate fluxes to go cc -> ec
  fluxes_pred(step, flux, x, B_cc);
  compute_E(step, Econv, x, B_cc, flux, 2);

  // add both electric fields together (glue the operators together)
  for (int p = 0; p < mrc_fld_nr_patches(E); p++) {
    mrc_fld_foreach(E, i,j,k, 2, 2) {
      for (int d=0; d < 3; d++) {
	M3(E, d, i,j,k, p) = M3(Ediff, d, i,j,k, p) + M3(Econv, d, i,j,k, p);
      }
    }  mrc_fld_foreach_end;
  }

  ggcm_mhd_put_3d_fld(mhd, B_cc);
  ggcm_mhd_put_3d_fld(mhd, Ediff);
  ggcm_mhd_put_3d_fld(mhd, Econv);
  ggcm_mhd_put_3d_fld(mhd, flux[0]);
  ggcm_mhd_put_3d_fld(mhd, flux[1]);
  ggcm_mhd_put_3d_fld(mhd, flux[2]);
  mrc_fld_put_as(E, Eout);
  // FIXME, should use _put_as, but don't want copy-back
  if (strcmp(mrc_fld_type(state_vec), FLD_TYPE) != 0) {
    mrc_fld_destroy(x);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_vlct, x)
static struct param ggcm_mhd_step_vlct_descr[] = {
  { "eqn"                , VAR(opt.eqn)            , PARAM_SELECT(OPT_EQN,
								  opt_eqn_descr)                },
  { "limiter"            , VAR(opt.limiter)        , PARAM_SELECT(OPT_LIMITER_GMINMOD,
								  opt_limiter_descr)            },
  { "riemann"            , VAR(opt.riemann)        , PARAM_SELECT(OPT_RIEMANN_HLL,
								  opt_riemann_descr)            },
  { "get_dt"             , VAR(opt.get_dt)         , PARAM_SELECT(OPT_GET_DT,
								  opt_get_dt_descr)             },
  { "background"         , VAR(opt.background)     , PARAM_BOOL(false)                          },

  { "debug_dump"         , VAR(debug_dump)         , PARAM_BOOL(false)                          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_ops

struct ggcm_mhd_step_ops ggcm_mhd_step_vlct_ops = {
  .name             = "vlct",
  .size             = sizeof(struct ggcm_mhd_step_vlct),
  .param_descr      = ggcm_mhd_step_vlct_descr,
  .setup            = ggcm_mhd_step_vlct_setup,
  .destroy          = ggcm_mhd_step_vlct_destroy,
  .run              = ggcm_mhd_step_vlct_run,
  .setup_flds       = ggcm_mhd_step_vlct_setup_flds,
  .get_e_ec         = ggcm_mhd_step_vlct_get_e_ec,
  .get_dt           = ggcm_mhd_step_vlct_get_dt,
};

