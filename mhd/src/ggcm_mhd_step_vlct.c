
#include <ggcm_mhd_step_private.h>

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <ggcm_mhd_diag.h>
#include <mhd_riemann.h>
#include <mhd_reconstruct.h>

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_io.h>

#include <mrc_fld_as_double.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include "mhd_1d.c"
#include "mhd_3d.c"


static const int nghost = 4; // FIXME

#include <string.h>

static int ldims[3];

// ======================================================================
// ggcm_mhd_step subclass "vlct"

struct ggcm_mhd_step_vlct {
  bool debug_dump;

  struct mhd_reconstruct *reconstruct_pred;
  struct mhd_reconstruct *reconstruct_corr;
  struct mhd_riemann *riemann;

  struct mrc_fld *U_1d;
  struct mrc_fld *U_l;
  struct mrc_fld *U_r;
  struct mrc_fld *W_1d;
  struct mrc_fld *W_l;
  struct mrc_fld *W_r;
  struct mrc_fld *F_1d;
  struct mrc_fld *Bxi;
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
  struct mrc_fld *Ecc = ggcm_mhd_step_get_3d_fld(step, 3);

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

  ggcm_mhd_step_put_3d_fld(step, Ecc);
}

static void
flux_pred(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc,
	  int ldim, int nghost, int j, int k, int dir, int p)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct mrc_fld *U_1d = sub->U_1d, *U_l = sub->U_l, *U_r = sub->U_r;
  struct mrc_fld *W_1d = sub->W_1d, *W_l = sub->W_l, *W_r = sub->W_r;
  struct mrc_fld *F_1d = sub->F_1d, *Bxi = sub->Bxi;

  pick_line_fc(U_1d, Bxi, x, B_cc, ldim, nghost, nghost, j, k, dir, p);
  mhd_prim_from_fc(step->mhd, W_1d, U_1d, ldim, nghost, nghost);
  mhd_reconstruct_run(sub->reconstruct_pred, U_l, U_r, W_l, W_r, W_1d, Bxi,
		      ldim, nghost - 1, nghost, dir);
  mhd_riemann_run(sub->riemann, F_1d, U_l, U_r, W_l, W_r, ldim, nghost - 1, nghost, dir);
  put_line_fc(flux[dir], F_1d, ldim, nghost - 1, nghost, j, k, dir, p);
}

static void
flux_corr(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc,
	  int ldim, int nghost, int j, int k, int dir, int p)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct mrc_fld *U_1d = sub->U_1d, *U_l = sub->U_l, *U_r = sub->U_r;
  struct mrc_fld *W_1d = sub->W_1d, *W_l = sub->W_l, *W_r = sub->W_r;
  struct mrc_fld *F_1d = sub->F_1d, *Bxi = sub->Bxi;

  pick_line_fc(U_1d, Bxi, x, B_cc, ldim, nghost - 1, nghost - 1, j, k, dir, p);
  mhd_prim_from_fc(step->mhd, W_1d, U_1d, ldim, nghost - 1, nghost - 1);
  mhd_reconstruct_run(sub->reconstruct_corr, U_l, U_r, W_l, W_r, W_1d, Bxi,
		      ldim, 99, 99, dir);
  mhd_riemann_run(sub->riemann, F_1d, U_l, U_r, W_l, W_r, ldim, 0, 1, dir);
  put_line_fc(flux[dir], F_1d, ldim, 0, 1, j, k, dir, p);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_create

static void
ggcm_mhd_step_vlct_create(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);

  mhd_reconstruct_set_type(sub->reconstruct_pred, "pcm_double");
  mhd_reconstruct_set_type(sub->reconstruct_corr, "plm_double");
  mhd_riemann_set_type(sub->riemann, "hll");
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_setup

static void
ggcm_mhd_step_vlct_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct ggcm_mhd *mhd = step->mhd;

  mhd_reconstruct_set_param_obj(sub->reconstruct_pred, "mhd", mhd);
  mhd_reconstruct_set_param_obj(sub->reconstruct_corr, "mhd", mhd);
  mhd_riemann_set_param_obj(sub->riemann, "mhd", mhd);

  setup_mrc_fld_1d(sub->U_1d, mhd->fld, 8);
  setup_mrc_fld_1d(sub->U_l , mhd->fld, 8);
  setup_mrc_fld_1d(sub->U_r , mhd->fld, 8);
  setup_mrc_fld_1d(sub->W_1d, mhd->fld, 8);
  setup_mrc_fld_1d(sub->W_l , mhd->fld, 8);
  setup_mrc_fld_1d(sub->W_r , mhd->fld, 8);
  setup_mrc_fld_1d(sub->F_1d, mhd->fld, 8);
  setup_mrc_fld_1d(sub->Bxi, mhd->fld, 1);

  assert(mhd->fld->_nr_ghosts >= 4);

  ggcm_mhd_step_setup_member_objs_sub(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_destroy

static void
ggcm_mhd_step_vlct_destroy(struct ggcm_mhd_step *step)
{
}

// ----------------------------------------------------------------------
// compute_Ediffu_const

static void
compute_Ediffu_const(struct ggcm_mhd_step *step, struct mrc_fld *E_ec,
                     struct mrc_fld *x, struct mrc_fld *B_cc)
{
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *j_ec = ggcm_mhd_step_get_3d_fld(step, 3);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  // calc edge-centered J
  for (int p = 0; p < mrc_fld_nr_patches(j_ec); p++) {
    double dx[3]; mrc_crds_get_dx(crds, p, dx);
    mrc_fld_data_t dxi[3] = { 1. / dx[0], 1. / dx[1], 1. / dx[2] };

    mrc_fld_foreach(j_ec, i,j,k, 0, 1) {
      M3(j_ec, 0, i,j,k, p) =
	(BZ_(x, i,j,k, p) - BZ_(x, i,j-1,k, p)) * dxi[1] -
	(BY_(x, i,j,k, p) - BY_(x, i,j,k-1, p)) * dxi[2];
      M3(j_ec, 1, i,j,k, p) =
	(BX_(x, i,j,k, p) - BX_(x, i,j,k-1, p)) * dxi[2] -
	(BZ_(x, i,j,k, p) - BZ_(x, i-1,j,k, p)) * dxi[0];
      M3(j_ec, 2, i,j,k, p) =
	(BY_(x, i,j,k, p) - BY_(x, i-1,j,k, p)) * dxi[0] -
	(BX_(x, i,j,k, p) - BX_(x, i,j-1,k, p)) * dxi[1];
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t eta = mhd->par.diffco / mhd->par.resnorm;

  for (int p = 0; p < mrc_fld_nr_patches(E_ec); p++) {
    mrc_fld_foreach(E_ec, i,j,k, 0, 1) {
      M3(E_ec, 0, i,j,k, p) = eta * M3(j_ec, 0, i,j,k, p);
      M3(E_ec, 1, i,j,k, p) = eta * M3(j_ec, 1, i,j,k, p);
      M3(E_ec, 2, i,j,k, p) = eta * M3(j_ec, 2, i,j,k, p);
    } mrc_fld_foreach_end;
  }

  // FIXME!!! energy flux

  ggcm_mhd_step_put_3d_fld(step, j_ec);
}

// ----------------------------------------------------------------------
// newstep_fc

// FIXME, take into account resistivity
// FIXME, rework

static mrc_fld_data_t
newstep_fc(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *Bcc)
{
  mrc_fld_data_t qsq,asq,cf1sq,cf2sq,cf3sq;
  mrc_fld_data_t b1,b2,b3,bsq,tsum,tdif;
  mrc_fld_data_t max_v1=0.0,max_v2=0.0,max_v3=0.0,max_dti = 0.0;

  mrc_fld_data_t gamma = mhd->par.gamm;
  mrc_fld_data_t gamma_minus_1 = gamma - 1.;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    double dx[3]; mrc_crds_get_dx(crds, p, dx);

    mrc_fld_foreach(x, i,j,k, 0, 0) {
      mrc_fld_data_t rri = 1.f / RR_(x, i,j,k, p);
      mrc_fld_data_t vx = RVX_(x, i,j,k, p) * rri;
      mrc_fld_data_t vy = RVY_(x, i,j,k, p) * rri;
      mrc_fld_data_t vz = RVZ_(x, i,j,k, p) * rri;
      qsq = sqr(vx) + sqr(vy) + sqr(vz);
      
      /* Use maximum of face-centered fields (always larger than cell-centered B) */
      b1 = M3(Bcc, 0, i,j,k, p) + fabs(BX_(x, i,j,k, p) - M3(Bcc, 0, i,j,k, p));
      b2 = M3(Bcc, 1, i,j,k, p) + fabs(BY_(x, i,j,k, p) - M3(Bcc, 1, i,j,k, p));
      b3 = M3(Bcc, 2, i,j,k, p) + fabs(BZ_(x, i,j,k, p) - M3(Bcc, 2, i,j,k, p));
      bsq = sqr(b1) + sqr(b2) + sqr(b3);
      /* compute sound speed squared */
      mrc_fld_data_t pp = fmax(gamma_minus_1*(EE_(x, i,j,k, p) - .5f*RR_(x, i,j,k, p)*qsq
					      - .5f*bsq), TINY_NUMBER);
      asq = gamma * pp * rri;
      
      /* compute fast magnetosonic speed squared in each direction */
      tsum = bsq * rri + asq;
      tdif = bsq * rri - asq;
      cf1sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b2*b2+b3*b3) * rri));
      cf2sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b3*b3) * rri));
      cf3sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b2*b2) * rri));
      
      /* compute maximum inverse dt (corresponding to minimum dt) */
      max_v1 = fmax(max_v1, (fabs(vx) + sqrt(cf1sq)) / dx[0]);
      max_v2 = fmax(max_v2, (fabs(vy) + sqrt(cf2sq)) / dx[1]);
      max_v3 = fmax(max_v3, (fabs(vz) + sqrt(cf3sq)) / dx[2]);
    } mrc_fld_foreach_end;
  }

  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);

  if (gdims[0] > 1) max_dti = fmax(max_dti, max_v1);
  if (gdims[1] > 1) max_dti = fmax(max_dti, max_v2);
  if (gdims[2] > 1) max_dti = fmax(max_dti, max_v3);

  mrc_fld_data_t cfl = mhd->par.thx;
  mrc_fld_data_t local_dt = cfl / max_dti;
  mrc_fld_data_t global_dt;
  MPI_Allreduce(&local_dt, &global_dt, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));
  return global_dt;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_run

static void
ggcm_mhd_step_vlct_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *B_cc = ggcm_mhd_step_get_3d_fld(step, 3);

  ldims[0] = mrc_fld_spatial_dims(x)[0];
  ldims[1] = mrc_fld_spatial_dims(x)[1];
  ldims[2] = mrc_fld_spatial_dims(x)[2];

  struct mrc_fld *x_half = ggcm_mhd_step_get_3d_fld(step, 8);
  struct mrc_fld *E_ec = ggcm_mhd_step_get_3d_fld(step, 3);
  struct mrc_fld *flux[3] = { ggcm_mhd_step_get_3d_fld(step, 8),
			      ggcm_mhd_step_get_3d_fld(step, 8),
			      ggcm_mhd_step_get_3d_fld(step, 8), };

  // CFL CONDITION

  ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);

  compute_B_cc(B_cc, x, 3, 3);
  if (step->do_nwst) {
    double old_dt = mhd->dt;
    mhd->dt = newstep_fc(mhd, x, B_cc);
    if (mhd->dt != old_dt) {
      mpi_printf(ggcm_mhd_comm(mhd), "switched dt %g <- %g\n", mhd->dt, old_dt);
      if (mhd->dt < mhd->par.dtmin) {
        mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin, aborting now!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }

  mrc_fld_data_t dt = mhd->dt;

  // resistivity

  if (mhd->par.magdiffu == MAGDIFFU_CONST) {
    if (mhd->par.diffco > 0.) {
      compute_Ediffu_const(step, E_ec, x, B_cc);
      update_ct_uniform(mhd, x, E_ec, dt, 0, 0);
      ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
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
    ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
    ggcm_mhd_diag_run_now(diag, x, DIAG_TYPE_3D, cnt++);
  }

  mrc_fld_copy_range(x_half, x, 0, 8);
  mhd_fluxes(step, flux, x, B_cc, nghost, nghost, flux_pred);
  update_finite_volume_uniform(mhd, x_half, flux, .5 * dt, nghost - 1, nghost - 1);
  compute_E(step, E_ec, x, B_cc, flux, 4);
  update_ct_uniform(mhd, x_half, E_ec, .5 * dt, nghost - 1, nghost - 1);

  // CORRECTOR

  compute_B_cc(B_cc, x_half, nghost - 1, nghost - 1);
  mhd_fluxes(step, flux, x_half, B_cc, 1, nghost, flux_corr);
  update_finite_volume_uniform(mhd, x, flux, dt, 0, 0);
  compute_E(step, E_ec, x_half, B_cc, flux, 4);
  update_ct_uniform(mhd, x, E_ec, dt, 0, 0);

  // clean up

  ggcm_mhd_step_put_3d_fld(step, B_cc);
  ggcm_mhd_step_put_3d_fld(step, x_half);
  ggcm_mhd_step_put_3d_fld(step, E_ec);
  ggcm_mhd_step_put_3d_fld(step, flux[0]);
  ggcm_mhd_step_put_3d_fld(step, flux[1]);
  ggcm_mhd_step_put_3d_fld(step, flux[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_setup_flds

static void
ggcm_mhd_step_vlct_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 4);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_FULLY_CONSERVATIVE);
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

  struct mrc_fld *B_cc = ggcm_mhd_step_get_3d_fld(step, 3);
  struct mrc_fld *Ediff = ggcm_mhd_step_get_3d_fld(step, 3);
  struct mrc_fld *Econv = ggcm_mhd_step_get_3d_fld(step, 3);
  struct mrc_fld *flux[3] = { ggcm_mhd_step_get_3d_fld(step, 8),
                              ggcm_mhd_step_get_3d_fld(step, 8),
                              ggcm_mhd_step_get_3d_fld(step, 8), };

  // Do diffusive terms first
  compute_B_cc(B_cc, x, 3, 3);
  if (mhd->par.magdiffu == MAGDIFFU_CONST) {
    compute_Ediffu_const(step, Ediff, x, B_cc);
  } else {
    mrc_fld_foreach(Ediff, i,j,k, 2, 2) {
    for (int d=0; d < 3; d++) {
        F3(Ediff, d, i,j,k) = 0.0;
      }
    }  mrc_fld_foreach_end;
    mpi_printf(ggcm_mhd_comm(mhd), "WARNING: magdiffu '%d' not handled for diags!\n",
               mhd->par.magdiffu);
  }

  // Do convective term using the apprepriate fluxes to go cc -> ec
  mhd_fluxes(step, flux, x, B_cc, nghost, nghost, flux_pred);
  compute_E(step, Econv, x, B_cc, flux, 2);

  // add both electric fields together (glue the operators together)
  mrc_fld_foreach(E, i,j,k, 2, 2) {
    for (int d=0; d < 3; d++) {
      F3(E, d, i,j,k) = F3(Ediff, d, i,j,k) + F3(Econv, d, i,j,k);
    }
  }  mrc_fld_foreach_end;

  ggcm_mhd_step_put_3d_fld(step, B_cc);
  ggcm_mhd_step_put_3d_fld(step, Ediff);
  ggcm_mhd_step_put_3d_fld(step, Econv);
  ggcm_mhd_step_put_3d_fld(step, flux[0]);
  ggcm_mhd_step_put_3d_fld(step, flux[1]);
  ggcm_mhd_step_put_3d_fld(step, flux[2]);
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
  { "debug_dump"      , VAR(debug_dump)      , PARAM_BOOL(false)            },

  { "reconstruct_pred", VAR(reconstruct_pred), MRC_VAR_OBJ(mhd_reconstruct) },
  { "reconstruct_corr", VAR(reconstruct_corr), MRC_VAR_OBJ(mhd_reconstruct) },
  { "riemann"         , VAR(riemann)         , MRC_VAR_OBJ(mhd_riemann)     },

  { "U_1d"            , VAR(U_1d)            , MRC_VAR_OBJ(mrc_fld)         },
  { "U_l"             , VAR(U_l)             , MRC_VAR_OBJ(mrc_fld)         },
  { "U_r"             , VAR(U_r)             , MRC_VAR_OBJ(mrc_fld)         },
  { "W_1d"            , VAR(W_1d)            , MRC_VAR_OBJ(mrc_fld)         },
  { "W_l"             , VAR(W_l)             , MRC_VAR_OBJ(mrc_fld)         },
  { "W_r"             , VAR(W_r)             , MRC_VAR_OBJ(mrc_fld)         },
  { "F_1d"            , VAR(F_1d)            , MRC_VAR_OBJ(mrc_fld)         },
  { "Bxi"             , VAR(Bxi)             , MRC_VAR_OBJ(mrc_fld)         },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct_ops

struct ggcm_mhd_step_ops ggcm_mhd_step_vlct_ops = {
  .name             = "vlct",
  .size             = sizeof(struct ggcm_mhd_step_vlct),
  .param_descr      = ggcm_mhd_step_vlct_descr,
  .create           = ggcm_mhd_step_vlct_create,
  .setup            = ggcm_mhd_step_vlct_setup,
  .destroy          = ggcm_mhd_step_vlct_destroy,
  .run              = ggcm_mhd_step_vlct_run,
  .setup_flds       = ggcm_mhd_step_vlct_setup_flds,
  .get_e_ec         = ggcm_mhd_step_vlct_get_e_ec,
};

