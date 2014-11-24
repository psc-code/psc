
#include <ggcm_mhd_step_private.h>

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <mhd_riemann.h>
#include <mhd_reconstruct.h>

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_fld_as_double_aos.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include "mhd_1d.c"
#include "mhd_3d.c"


static const int nghost = 4; // FIXME

static int ldims[3];

// ======================================================================
// ggcm_mhd_step subclass "vlct"

struct ggcm_mhd_step_vlct {
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

#define MAKE_COMPUTE_E_EDGE(X,Y,Z,IJ,JJ,KJ,IK,JK,KK)			\
									\
static void								\
compute_E##X##_edge(struct mrc_fld *E, struct mrc_fld *Ecc,		\
		    struct mrc_fld **fluxes, int bnd)			\
{									\
  int BY = BX + Y, BZ = BX + Z;						\
									\
  mrc_fld_foreach(E, i,j,k, bnd, bnd+1) {				\
    mrc_fld_data_t de1_l2, de1_r2, de1_l3, de1_r3;			\
									\
    de1_l3 = l_r_avg(F3(fluxes[Y], RR, i-IK,j-JK,k-KK),			\
		     F3(fluxes[Z], BY, i   ,j   ,k   ) - F3(Ecc, X, i   -IK,j   -JK,k   -KK), \
		     F3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ) - F3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK)); \
    									\
    de1_r3 = l_r_avg(F3(fluxes[Y], RR, i,j,k),				\
		     F3(fluxes[Z], BY, i   ,j   ,k   ) - F3(Ecc, X, i   ,j   ,k   ), \
		     F3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ) - F3(Ecc, X, i-IJ,j-JJ,k-KJ)); \
    									\
    de1_l2 = l_r_avg( F3(fluxes[Z], RR, i-IJ,j-JJ,k-KJ),			\
		     -F3(fluxes[Y], BZ, i   ,j   ,k   ) - F3(Ecc, X, i-IJ   ,j-JJ   ,k-KJ   ), \
		     -F3(fluxes[Y], BZ, i-IK,j-JK,k-KK) - F3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK)); \
    									\
    de1_r2 = l_r_avg( F3(fluxes[Z], RR, i,j,k),				\
		     -F3(fluxes[Y], BZ, i   ,j   ,k   ) - F3(Ecc, X, i   ,j   ,k   ), \
		     -F3(fluxes[Y], BZ, i-IK,j-JK,k-KK) - F3(Ecc, X, i-IK,j-JK,k-KK)); \
    									\
    F3(E, X, i,j,k) =						\
      .25 * (  F3(fluxes[Z], BY, i,j,k) + F3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ) \
	     - F3(fluxes[Y], BZ, i,j,k) - F3(fluxes[Y], BZ, i-IK,j-JK,k-KK) \
	     + de1_l2 + de1_r2 + de1_l3 + de1_r3);		        \
  } mrc_fld_foreach_end;						\
}

#define MAKE_COMPUTE_E_EDGE_(X,Y,Z,IJ,JJ,KJ,IK,JK,KK)			\
									\
static void								\
compute_E##X##_edge(struct mrc_fld *E, struct mrc_fld *Ecc,		\
		   struct mrc_fld **fluxes, int bnd)			\
{									\
  int BY = BX + Y, BZ = BX + Z;						\
									\
  mrc_fld_foreach(E, i,j,k, bnd, bnd+1) {				\
    mrc_fld_data_t de1_l2, de1_r2, de1_l3, de1_r3;			\
									\
    de1_l3 = l_r_avg(F3(fluxes[Y], RR, i-IK,j-JK,k-KK),			\
		     F3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ) - F3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK), \
		     F3(fluxes[Z], BY, i   ,j   ,k   ) - F3(Ecc, X, i   -IK,j   -JK,k   -KK)); \
    									\
    de1_r3 = l_r_avg(F3(fluxes[Y], RR, i,j,k),				\
		     F3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ) - F3(Ecc, X, i-IJ,j-JJ,k-KJ), \
		     F3(fluxes[Z], BY, i   ,j   ,k   ) - F3(Ecc, X, i   ,j   ,k   )); \
    									\
    de1_l2 = l_r_avg( F3(fluxes[Z], RR, i-IJ,j-JJ,k-KJ),			\
		     -F3(fluxes[Y], BZ, i-IK,j-JK,k-KK) - F3(Ecc, X, i-IJ-IK,j-JJ-JK,k-KJ-KK), \
		     -F3(fluxes[Y], BZ, i   ,j   ,k   ) - F3(Ecc, X, i-IJ   ,j-JJ   ,k-KJ   )); \
    									\
    de1_r2 = l_r_avg( F3(fluxes[Z], RR, i,j,k),				\
		     -F3(fluxes[Y], BZ, i-IK,j-JK,k-KK) - F3(Ecc, X, i-IK,j-JK,k-KK), \
		     -F3(fluxes[Y], BZ, i   ,j   ,k   ) - F3(Ecc, X, i   ,j   ,k   )); \
    									\
    F3(E, X, i,j,k) =						\
      .25 * (  F3(fluxes[Z], BY, i,j,k) + F3(fluxes[Z], BY, i-IJ,j-JJ,k-KJ) \
	     - F3(fluxes[Y], BZ, i,j,k) - F3(fluxes[Y], BZ, i-IK,j-JK,k-KK) \
	     + de1_l2 + de1_r2 + de1_l3 + de1_r3);		        \
  } mrc_fld_foreach_end;						\
}

MAKE_COMPUTE_E_EDGE(0,1,2, 0,1,0, 0,0,1)
MAKE_COMPUTE_E_EDGE(1,2,0, 0,0,1, 1,0,0)
MAKE_COMPUTE_E_EDGE_(2,0,1, 1,0,0, 0,1,0)

static void
compute_E(struct ggcm_mhd_step *step, struct mrc_fld *E,
	  struct mrc_fld *x, struct mrc_fld *B_cc, struct mrc_fld *fluxes[3],
	  int bnd)
{
  struct mrc_fld *Ecc = ggcm_mhd_step_get_3d_fld(step, 3);

  // calculate cell-centered E first
  mrc_fld_foreach(Ecc, i,j,k, bnd, bnd) {
    mrc_fld_data_t rri = 1. / F3(x, RR, i,j,k);
    mrc_fld_data_t B[3] = { F3(B_cc, 0, i,j,k), F3(B_cc, 1, i,j,k), F3(B_cc, 2, i,j,k) };
    mrc_fld_data_t v[3] = { rri * F3(x, RVX, i,j,k), rri * F3(x, RVY, i,j,k), rri * F3(x, RVZ, i,j,k) };
    F3(Ecc, 0, i,j,k) = B[1] * v[2] - B[2] * v[1];
    F3(Ecc, 1, i,j,k) = B[2] * v[0] - B[0] * v[2];
    F3(Ecc, 2, i,j,k) = B[0] * v[1] - B[1] * v[0];
  } mrc_fld_foreach_end;

  // then calculate edge-centered E based on the cell-centered one and
  // the fluxes
  compute_E0_edge(E, Ecc, fluxes, bnd - 1);
  compute_E1_edge(E, Ecc, fluxes, bnd - 1);
  compute_E2_edge(E, Ecc, fluxes, bnd - 1);

  ggcm_mhd_step_put_3d_fld(step, Ecc);
}

static void
flux_pred(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc,
	  int ldim, int nghost, int j, int k, int dir)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct mrc_fld *U_1d = sub->U_1d, *U_l = sub->U_l, *U_r = sub->U_r;
  struct mrc_fld *W_1d = sub->W_1d, *W_l = sub->W_l, *W_r = sub->W_r;
  struct mrc_fld *F_1d = sub->F_1d, *Bxi = sub->Bxi;

  pick_line_fc(U_1d, Bxi, x, B_cc, ldim, nghost, nghost, j, k, dir);
  mhd_prim_from_fc(step->mhd, W_1d, U_1d, ldim, nghost, nghost);
  mhd_reconstruct_run(sub->reconstruct_pred, U_l, U_r, W_l, W_r, W_1d, Bxi,
		      ldim, nghost - 1, nghost, dir);
  mhd_riemann_run(sub->riemann, F_1d, U_l, U_r, W_l, W_r, ldim, nghost - 1, nghost, dir);
  put_line_fc(flux[dir], F_1d, ldim, nghost - 1, nghost, j, k, dir);
}

static void
flux_corr(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc,
	  int ldim, int nghost, int j, int k, int dir)
{
  struct ggcm_mhd_step_vlct *sub = ggcm_mhd_step_vlct(step);
  struct mrc_fld *U_1d = sub->U_1d, *U_l = sub->U_l, *U_r = sub->U_r;
  struct mrc_fld *W_1d = sub->W_1d, *W_l = sub->W_l, *W_r = sub->W_r;
  struct mrc_fld *F_1d = sub->F_1d, *Bxi = sub->Bxi;

  pick_line_fc(U_1d, Bxi, x, B_cc, ldim, nghost - 1, nghost - 1, j, k, dir);
  mhd_prim_from_fc(step->mhd, W_1d, U_1d, ldim, nghost - 1, nghost - 1);
  mhd_reconstruct_run(sub->reconstruct_corr, U_l, U_r, W_l, W_r, W_1d, Bxi,
		      ldim, 99, 99, dir);
  mhd_riemann_run(sub->riemann, F_1d, U_l, U_r, W_l, W_r, ldim, 0, 1, dir);
  put_line_fc(flux[dir], F_1d, ldim, 0, 1, j, k, dir);
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
// magdiffu_const

static void
magdiffu_const(struct ggcm_mhd_step *step, struct mrc_fld *x, struct mrc_fld *B_cc,
	       mrc_fld_data_t dt)
{
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *j_ec = ggcm_mhd_step_get_3d_fld(step, 3);
  struct mrc_fld *E_ec = ggcm_mhd_step_get_3d_fld(step, 3);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double dx[3]; mrc_crds_get_dx(crds, dx);
  mrc_fld_data_t dxi[3] = { 1. / dx[0], 1. / dx[1], 1. / dx[2] };

  // calc edge-centered J
  mrc_fld_foreach(j_ec, i,j,k, 0, 1) {
    F3(j_ec, 0, i,j,k) =
      (BZ(x, i,j,k) - BZ(x, i,j-1,k)) * dxi[1] -
      (BY(x, i,j,k) - BY(x, i,j,k-1)) * dxi[2];
    F3(j_ec, 1, i,j,k) =
      (BX(x, i,j,k) - BX(x, i,j,k-1)) * dxi[2] -
      (BZ(x, i,j,k) - BZ(x, i-1,j,k)) * dxi[0];
    F3(j_ec, 2, i,j,k) =
      (BY(x, i,j,k) - BY(x, i-1,j,k)) * dxi[0] -
      (BX(x, i,j,k) - BX(x, i,j-1,k)) * dxi[1];
  } mrc_fld_foreach_end;

  mrc_fld_data_t eta = mhd->par.diffco / mhd->par.resnorm;

  mrc_fld_foreach(j_ec, i,j,k, 0, 1) {
    F3(E_ec, 0, i,j,k) = eta * F3(j_ec, 0, i,j,k);
    F3(E_ec, 1, i,j,k) = eta * F3(j_ec, 1, i,j,k);
    F3(E_ec, 2, i,j,k) = eta * F3(j_ec, 2, i,j,k);
  } mrc_fld_foreach_end;

  // FIXME!!! energy flux

  update_ct_uniform(mhd, x, E_ec, dt, 0, 0);

  ggcm_mhd_step_put_3d_fld(step, j_ec);
  ggcm_mhd_step_put_3d_fld(step, E_ec);
}

// ----------------------------------------------------------------------
// newstep_fc

// FIXME, take into account resistivity

static mrc_fld_data_t
newstep_fc(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *Bcc)
{
  mrc_fld_data_t di,v1,v2,v3,qsq,asq,cf1sq,cf2sq,cf3sq;
  mrc_fld_data_t p;
  mrc_fld_data_t b1,b2,b3,bsq,tsum,tdif;
  mrc_fld_data_t max_v1=0.0,max_v2=0.0,max_v3=0.0,max_dti = 0.0;

  mrc_fld_data_t gamma = mhd->par.gamm;
  mrc_fld_data_t gamma_minus_1 = gamma - 1.;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  double dx[3]; mrc_crds_get_dx(crds, dx);

  mrc_fld_foreach(x, i,j,k, 0, 0) {
    di = 1.0/RR(x, i,j,k);
    v1 = RVX(x, i,j,k)*di;
    v2 = RVY(x, i,j,k)*di;
    v3 = RVZ(x, i,j,k)*di;
    qsq = v1*v1 + v2*v2 + v3*v3;
    
    /* Use maximum of face-centered fields (always larger than cell-centered B) */
    b1 = F3(Bcc, 0, i,j,k) + fabs(BX(x, i,j,k) - F3(Bcc, 0, i,j,k));
    b2 = F3(Bcc, 1, i,j,k) + fabs(BY(x, i,j,k) - F3(Bcc, 1, i,j,k));
    b3 = F3(Bcc, 2, i,j,k) + fabs(BZ(x, i,j,k) - F3(Bcc, 2, i,j,k));
    bsq = sqr(b1) + sqr(b2) + sqr(b3);
    /* compute sound speed squared */
    p = fmax(gamma_minus_1*(EE(x, i,j,k) - 0.5*RR(x, i,j,k)*qsq
		      - 0.5*bsq), TINY_NUMBER);
    asq = gamma*p*di;
    
    /* compute fast magnetosonic speed squared in each direction */
    tsum = bsq*di + asq;
    tdif = bsq*di - asq;
    cf1sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b2*b2+b3*b3)*di));
    cf2sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b3*b3)*di));
    cf3sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b2*b2)*di));
    
    /* compute maximum cfl velocity (corresponding to minimum dt) */
    max_v1 = fmax(max_v1, fabs(v1) + sqrt(cf1sq));
    max_v2 = fmax(max_v2, fabs(v2) + sqrt(cf2sq));
    max_v3 = fmax(max_v3, fabs(v3) + sqrt(cf3sq));
    
  } mrc_fld_foreach_end;
  
  max_dti = fmax(max_dti, max_v1/dx[0]);
  max_dti = fmax(max_dti, max_v2/dx[1]);
  max_dti = fmax(max_dti, max_v3/dx[2]);

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
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *B_cc = ggcm_mhd_step_get_3d_fld(step, 3);

  ldims[0] = mrc_fld_spatial_dims(x)[0];
  ldims[1] = mrc_fld_spatial_dims(x)[1];
  ldims[2] = mrc_fld_spatial_dims(x)[2];

  struct mrc_fld *x_half = ggcm_mhd_step_get_3d_fld(step, 8);
  struct mrc_fld *E = ggcm_mhd_step_get_3d_fld(step, 3);
  struct mrc_fld *flux[3] = { ggcm_mhd_step_get_3d_fld(step, 8),
			      ggcm_mhd_step_get_3d_fld(step, 8),
			      ggcm_mhd_step_get_3d_fld(step, 8), };

  // CFL CONDITION

  ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);

  compute_B_cc(B_cc, x, 3, 3);
  if (step->do_nwst) {
    double old_dt = mhd->dt;
    mhd->dt = newstep_fc(mhd, x, B_cc);
    mpi_printf(ggcm_mhd_comm(mhd), "switched dt %g <- %g\n", mhd->dt, old_dt);
    if (mhd->dt < mhd->par.dtmin) {  
      mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin, aborting now!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  mrc_fld_data_t dt = mhd->dt;

  // resistivity

  if (mhd->par.magdiffu == MAGDIFFU_CONST) {
    magdiffu_const(step, x, B_cc, dt);
    ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
    compute_B_cc(B_cc, x, 3, 3);
  } else {
    mpi_printf(ggcm_mhd_comm(mhd), "WARNING: magdiffu '%d' not handled!\n",
	       mhd->par.magdiffu);
  }

  // PREDICTOR

  mrc_fld_copy_range(x_half, x, 0, 8);
  mhd_fluxes(step, flux, x, B_cc, nghost, nghost, flux_pred);
  update_finite_volume_uniform(mhd, x_half, flux, .5 * dt, nghost - 1, nghost - 1);
  compute_E(step, E, x, B_cc, flux, 4);
  update_ct_uniform(mhd, x_half, E, .5 * dt, nghost - 1, nghost - 1);

  // CORRECTOR

  compute_B_cc(B_cc, x_half, nghost - 1, nghost - 1);
  mhd_fluxes(step, flux, x_half, B_cc, 1, nghost, flux_corr);
  update_finite_volume_uniform(mhd, x, flux, dt, 0, 0);
  compute_E(step, E, x_half, B_cc, flux, 4);
  update_ct_uniform(mhd, x, E, dt, 0, 0);

  // clean up

  ggcm_mhd_step_put_3d_fld(step, B_cc);
  ggcm_mhd_step_put_3d_fld(step, x_half);
  ggcm_mhd_step_put_3d_fld(step, E);
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
  mrc_fld_set_param_int(mhd->fld, "nr_comps", _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vlct subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_vlct, x)
static struct param ggcm_mhd_step_vlct_descr[] = {
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
};

