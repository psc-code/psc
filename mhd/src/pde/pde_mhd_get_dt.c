
#include "ggcm_mhd_crds.h" // FIXME

#ifdef FD1X

// ----------------------------------------------------------------------
// patch_get_dt_scons_ggcm_c

static mrc_fld_data_t
patch_get_dt_scons_ggcm_c(fld3d_t p_f)
{
  fld3d_t p_W, p_U, p_cmsv, p_bcc, p_ymask, p_zmask;
  fld3d_setup_view(&p_W    , p_f, _RR);
  fld3d_setup_view(&p_U    , p_f, _RR1);
  fld3d_setup_view(&p_cmsv , p_f, _CMSV);
  fld3d_setup_view(&p_bcc  , p_f, _BX);
  fld3d_setup_view(&p_ymask, p_f, _YMASK);
  fld3d_setup_view(&p_zmask, p_f, _ZMASK);
  
  patch_primvar(p_W, p_U, p_cmsv);
  patch_primbb(p_bcc, p_U);
  patch_zmaskn(p_zmask, p_W, p_bcc, p_ymask);

  mrc_fld_data_t dt = 1e10f;

  mrc_fld_data_t splim = s_speedlimit_code;
  mrc_fld_data_t eps   = 1e-9f;

  fld3d_foreach(i, j, k, 0, 0) {
    mrc_fld_data_t hh = mrc_fld_max(mrc_fld_max(FD1X(i), FD1Y(j)), FD1Z(k));
    mrc_fld_data_t rri = 1.f / mrc_fld_abs(F3S(p_f, RR, i,j,k)); // FIME abs necessary?
    mrc_fld_data_t bb = sqr(F3S(p_f, _BX, i,j,k)) + sqr(F3S(p_f, _BY, i,j,k)) + sqr(F3S(p_f, _BZ, i,j,k));
    mrc_fld_data_t vv1 = fminf(sqrtf(bb * rri), splim);
    mrc_fld_data_t pp = F3S(p_f, _PP, i,j,k);
    mrc_fld_data_t vv2 = sqrtf(s_gamma * mrc_fld_max(0.f, pp) * rri);
    mrc_fld_data_t vv3 = sqrtf(sqr(F3S(p_f, _VX, i,j,k)) + sqr(F3S(p_f, _VY, i,j,k)) + sqr(F3S(p_f, _VZ, i,j,k)));
    mrc_fld_data_t vv = sqrtf(sqr(vv1) + sqr(vv2)) + vv3;
    vv = mrc_fld_max(eps, vv);
    
    mrc_fld_data_t zmask = F3S(p_f, _ZMASK, i,j,k);
    mrc_fld_data_t tt = s_cfl / mrc_fld_max(eps, hh*vv*zmask);
    dt = mrc_fld_min(dt, tt);
  } fld3d_foreach_end;

  return dt;
}

// ----------------------------------------------------------------------
// patch_get_dt_scons_ggcm_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define newstep_F77 F77_FUNC(newstep,NEWSTEP)

void newstep_F77(real *pp, real *rr, real *vx, real *vy, real *vz,
		 real *bx, real *by, real *bz, real *zmask, real *dtn);

static mrc_fld_data_t
patch_get_dt_scons_ggcm_fortran(fld3d_t p_f)
{
  fld3d_t p_W, p_U, p_cmsv, p_bcc, p_ymask, p_zmask;
  fld3d_setup_view(&p_W    , p_f, _RR);
  fld3d_setup_view(&p_U    , p_f, _RR1);
  fld3d_setup_view(&p_cmsv , p_f, _CMSV);
  fld3d_setup_view(&p_bcc  , p_f, _BX);
  fld3d_setup_view(&p_ymask, p_f, _YMASK);
  fld3d_setup_view(&p_zmask, p_f, _ZMASK);
  
  patch_primvar(p_W, p_U, p_cmsv);
  patch_primbb(p_bcc, p_U);
  patch_zmaskn(p_zmask, p_W, p_bcc, p_ymask);

  real dtn;
  newstep_F77(F(p_f, _PP), F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ),
	      F(p_f, _BX), F(p_f, _BY), F(p_f, _BZ), F(p_f, _ZMASK), &dtn);

  return dtn;
}

#endif // HAVE_OPENGGCM && MRC_FLD_AS_FLOAT_H

// ----------------------------------------------------------------------
// patch_get_dt_scons_ggcm

static mrc_fld_data_t
patch_get_dt_scons_ggcm(fld3d_t p_f)
{
  if (s_opt_mhd_newstep == OPT_MHD_C) {
    return patch_get_dt_scons_ggcm_c(p_f);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_newstep == OPT_MHD_FORTRAN) {
    return patch_get_dt_scons_ggcm_fortran(p_f);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// pde_mhd_get_dt_scons_ggcm

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_scons_ggcm(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  fld3d_t p_f;
  fld3d_setup(&p_f, x);

  mrc_fld_data_t dt = 1e10f;
  pde_for_each_patch(p) {
    fld3d_get(&p_f, p);
    dt = mrc_fld_min(dt, patch_get_dt_scons_ggcm(p_f));
    fld3d_put(&p_f, p);
  }
  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));
  
  return dtn;
}

#endif

// ----------------------------------------------------------------------
// pde_mhd_get_dt_scons

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_scons(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *zmask, 
		     int m_zmask)
{
  assert(s_opt_eqn == OPT_EQN_MHD_SCONS);

  struct mrc_fld *b0 = mhd->b0;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  mrc_fld_data_t splim2 = sqr(mhd->par.speedlimit / mhd->vvnorm);
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;
  mrc_fld_data_t d_i    = mhd->par.d_i;

  mrc_fld_data_t thx    = mhd->par.thx;
  mrc_fld_data_t eps    = 1e-9f;
  mrc_fld_data_t dt     = 1e10f;

  mrc_fld_data_t two_pi_d_i = 2. * M_PI * d_i;
  bool have_hall = d_i > 0.f;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    float *fd1x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, FD1, p);
    float *fd1y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, FD1, p);
    float *fd1z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, FD1, p);

    mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
      mrc_fld_data_t hh = mrc_fld_max(mrc_fld_max(fd1x[ix], fd1y[iy]), fd1z[iz]);
      mrc_fld_data_t rri = 1.f / mrc_fld_abs(RR_(x, ix,iy,iz, p)); // FIXME abs necessary?
      mrc_fld_data_t bb = (sqr(.5f * (BTX_(x, ix,iy,iz, p) + BTX_(x, ix+dx,iy,iz, p))) + 
			   sqr(.5f * (BTY_(x, ix,iy,iz, p) + BTY_(x, ix,iy+dy,iz, p))) +
			   sqr(.5f * (BTZ_(x, ix,iy,iz, p) + BTZ_(x, ix,iy,iz+dz, p))));
      mrc_fld_data_t rrvv = (sqr(RVX_(x, ix,iy,iz, p)) + 
			     sqr(RVY_(x, ix,iy,iz, p)) +
			     sqr(RVZ_(x, ix,iy,iz, p)));
      mrc_fld_data_t pp = gamma_m1 * (UU_(x, ix,iy,iz, p) - .5f * rrvv * rri);

      if (have_hall) {
	bb *= 1 + sqr(two_pi_d_i * hh);
      }      
      
      mrc_fld_data_t vA2 = mrc_fld_min(bb * rri, splim2);
      mrc_fld_data_t cs2 = s_gamma * pp * rri;
      mrc_fld_data_t vv = mrc_fld_sqrt(vA2 + cs2) + mrc_fld_sqrt(rrvv) * rri;
      vv = mrc_fld_max(eps, vv);
      
      mrc_fld_data_t zm = 1.f;
      if (zmask) {
	zm = M3(zmask, m_zmask, ix,iy,iz, p);
      }
      mrc_fld_data_t tt = thx / mrc_fld_max(eps, hh*vv*zm);
      dt = mrc_fld_min(dt, tt);
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  return dtn;
}

// ----------------------------------------------------------------------
// pde_mhd_get_dt_fcons

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_fcons(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  assert(!s_opt_background);

  static fld1d_state_t V, U;
  if (!fld1d_state_is_setup(V)) {
    fld1d_state_setup(&V);
    fld1d_state_setup(&U);
  }

  mrc_fld_data_t inv_dt = 0.f;
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);

    pde_for_each_dir(dir) {
      pde_line_set_dir(dir);
      int ib = 0, ie = s_ldims[dir];
      pde_for_each_line(dir, j, k, 0) {
	mhd_get_line_state(U, x, j, k, dir, p, ib, ie);
	mhd_prim_from_cons(V, U, ib, ie);
	for (int i = ib; i < ie; i++) {
	  mrc_fld_data_t *v = &F1S(V, 0, i);
	  // This is iffy: we should call wavespeed_mhd_fcons even if we're using
	  // scons variables, since we want all MHD waves taken into account, not just
	  // the sound waves. FIXME, there must be a better way
	  mrc_fld_data_t cf = wavespeed_mhd_fcons(NULL, v, i);

	  if (s_sw[dir]) inv_dt = mrc_fld_max(inv_dt, (mrc_fld_abs(v[VX]) + cf) * PDE_INV_DS(i));
	}
      }
    }
  }
  mrc_fld_data_t dt = mhd->par.thx / inv_dt;

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  if (s_opt_divb == OPT_DIVB_GLM) {
    s_divb_glm_ch = s_divb_glm_ch_fac * mhd->par.thx * s_g_dxyzmin / dtn;
  }

  return dtn;
}

#ifdef OPT_DIVB_CT

// ----------------------------------------------------------------------
// pde_mhd_get_dt_fcons_ct

// FIXME, take into account resistivity
// FIXME, rework

static void compute_B_cc(struct mrc_fld *B_cc, struct mrc_fld *x, int l, int r);


static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_fcons_ct(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  assert(s_opt_eqn == OPT_EQN_MHD_FCONS);

  static const mrc_fld_data_t eps = 1.e-10; // FIXME

  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;
  mrc_fld_data_t d_i = mhd->par.d_i;

  struct mrc_fld *Bcc = ggcm_mhd_get_3d_fld(mhd, 3);
  ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
  compute_B_cc(Bcc, x, 0, 0);

  mrc_fld_data_t max_dti_x = 0., max_dti_y = 0., max_dti_z = 0.;
  mrc_fld_data_t max_dti_diff = 0.;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);

    mrc_fld_foreach(x, i,j,k, 0, 0) {
      mrc_fld_data_t rri = 1.f / RR_(x, i,j,k, p);
      mrc_fld_data_t vx = RVX_(x, i,j,k, p) * rri;
      mrc_fld_data_t vy = RVY_(x, i,j,k, p) * rri;
      mrc_fld_data_t vz = RVZ_(x, i,j,k, p) * rri;
      mrc_fld_data_t vv = sqr(vx) + sqr(vy) + sqr(vz);
      
      /* Use maximum of face-centered fields (always larger than cell-centered B) */
      mrc_fld_data_t bx = M3(Bcc, 0, i,j,k, p) + fabs(BX_(x, i,j,k, p) - M3(Bcc, 0, i,j,k, p));
      mrc_fld_data_t by = M3(Bcc, 1, i,j,k, p) + fabs(BY_(x, i,j,k, p) - M3(Bcc, 1, i,j,k, p));
      mrc_fld_data_t bz = M3(Bcc, 2, i,j,k, p) + fabs(BZ_(x, i,j,k, p) - M3(Bcc, 2, i,j,k, p));
      mrc_fld_data_t bb = sqr(bx) + sqr(by) + sqr(bz);
      mrc_fld_data_t pp = mrc_fld_max(gamma_m1 * (EE_(x, i,j,k, p) - .5f*RR_(x, i,j,k, p)*vv
					      - .5f*bb), eps);
      mrc_fld_data_t cs2 = s_gamma * pp * rri;
      
      /* compute fast magnetosonic speed squared in each direction */
      mrc_fld_data_t tsum = bb * rri + cs2;
      mrc_fld_data_t tdif = bb * rri - cs2;
      mrc_fld_data_t cf1sq = .5f * (tsum + sqrt(sqr(tdif) + 4.f * cs2 * (sqr(by) + sqr(bz)) * rri));
      mrc_fld_data_t cf2sq = .5f * (tsum + sqrt(sqr(tdif) + 4.f * cs2 * (sqr(bx) + sqr(bz)) * rri));
      mrc_fld_data_t cf3sq = .5f * (tsum + sqrt(sqr(tdif) + 4.f * cs2 * (sqr(bx) + sqr(by)) * rri));
      
      max_dti_x = mrc_fld_max(max_dti_x, (mrc_fld_abs(vx) + mrc_fld_sqrt(cf1sq)) * PDE_INV_DX(i));
      max_dti_y = mrc_fld_max(max_dti_y, (mrc_fld_abs(vy) + mrc_fld_sqrt(cf2sq)) * PDE_INV_DY(j));
      max_dti_z = mrc_fld_max(max_dti_z, (mrc_fld_abs(vz) + mrc_fld_sqrt(cf3sq)) * PDE_INV_DZ(k));

      if (d_i > 0.) {
	mrc_fld_data_t inv_dx = mrc_fld_max(PDE_INV_DX(i), mrc_fld_max(PDE_INV_DY(j), PDE_INV_DZ(k)));
	
	// FIXME? in athena it was (dxmin**2 / 6), but i don't think that's right?
	max_dti_diff = mrc_fld_max(max_dti_diff, (d_i * bb / RR_(x, i,j,k, p)) * 16.f * sqr(inv_dx));
      }
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t max_dti = 0.;
  if (s_sw[0]) max_dti = mrc_fld_max(max_dti, max_dti_x);
  if (s_sw[1]) max_dti = mrc_fld_max(max_dti, max_dti_y);
  if (s_sw[2]) max_dti = mrc_fld_max(max_dti, max_dti_z);
  max_dti = mrc_fld_max(max_dti, max_dti_diff);

  mrc_fld_data_t dt = mhd->par.thx / max_dti;
  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  ggcm_mhd_put_3d_fld(mhd, Bcc);

  return dtn;
}

#endif

// ----------------------------------------------------------------------
// pde_mhd_get_dt_hd

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_hd(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  assert(s_opt_eqn == OPT_EQN_HD);

  static const mrc_fld_data_t eps = 1.e-10; // FIXME

  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  mrc_fld_data_t max_dti = 0.;
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    double dx[3]; mrc_crds_get_dx(crds, p, dx);

    mrc_fld_foreach(x, i,j,k, 0, 0) {
      mrc_fld_data_t rri = 1.f / RR_(x, i,j,k, p);
      mrc_fld_data_t vx = RVX_(x, i,j,k, p) * rri;
      mrc_fld_data_t vy = RVY_(x, i,j,k, p) * rri;
      mrc_fld_data_t vz = RVZ_(x, i,j,k, p) * rri;
      mrc_fld_data_t vv = sqr(vx) + sqr(vy) + sqr(vz);
      
      /* compute sound speed */
      mrc_fld_data_t pp = mrc_fld_max(gamma_m1 * (EE_(x, i,j,k, p) - .5f * RR_(x, i,j,k, p) * vv), eps);
      mrc_fld_data_t cs = mrc_fld_sqrt(s_gamma * pp * rri);
      
      /* compute min dt based on maximum wave velocity */
      if (s_sw[0]) { max_dti = mrc_fld_max(max_dti, (mrc_fld_abs(vx) + cs) / PDE_INV_DX(i)); }
      if (s_sw[1]) { max_dti = mrc_fld_max(max_dti, (mrc_fld_abs(vy) + cs) / PDE_INV_DY(j)); }
      if (s_sw[2]) { max_dti = mrc_fld_max(max_dti, (mrc_fld_abs(vz) + cs) / PDE_INV_DZ(k)); }
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t dt = mhd->par.thx / max_dti;
  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  return dtn;
}

// ----------------------------------------------------------------------
// pde_mhd_get_dt

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  if (s_opt_get_dt == OPT_GET_DT_MHD_GGCM) {
    assert(0);
  } else if (s_opt_get_dt == OPT_GET_DT_MHD) { // FIXME, name etc
    return pde_mhd_get_dt_fcons(mhd, x);
#ifdef OPT_DIVB_CT
  } else if (s_opt_get_dt == OPT_GET_DT_MHD_CT) {
    return pde_mhd_get_dt_fcons_ct(mhd, x);
#endif
  } else if (s_opt_get_dt == OPT_GET_DT_HD) {
    return pde_mhd_get_dt_hd(mhd, x);
  } else {
    assert(0);
  }
}
