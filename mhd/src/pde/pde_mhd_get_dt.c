
#include "ggcm_mhd_crds.h" // FIXME

// ----------------------------------------------------------------------
// pde_mhd_get_dt_scons

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_scons(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *zmask, 
		     int m_zmask)
{
  struct mrc_fld *b0 = mhd->b0;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] > 1), dy = (gdims[1] > 1), dz = (gdims[2] > 1);

  mrc_fld_data_t splim2 = sqr(mhd->par.speedlimit / mhd->vvnorm);
  mrc_fld_data_t gamm   = mhd->par.gamm;
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
      if (have_hall) {
	bb *= 1 + sqr(two_pi_d_i * hh);
      }      
      
      mrc_fld_data_t vv1 = mrc_fld_min(bb * rri, splim2);
      mrc_fld_data_t rrvv = (sqr(RVX_(x, ix,iy,iz, p)) + 
			     sqr(RVY_(x, ix,iy,iz, p)) +
			     sqr(RVZ_(x, ix,iy,iz, p)));
      mrc_fld_data_t pp = (gamm - 1.f) * (UU_(x, ix,iy,iz, p) - .5f * rrvv * rri);
      mrc_fld_data_t vv2 = gamm * pp * rri;
      mrc_fld_data_t vv = mrc_fld_sqrt(vv1 + vv2) + mrc_fld_sqrt(rrvv) * rri;
      vv = mrc_fld_max(eps, vv);
      
      mrc_fld_data_t zm = 1.f;
      if (zmask) {
	zm = F3(zmask, m_zmask, ix,iy,iz);
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
  struct mrc_fld *b0 = mhd->b0;

  mrc_fld_data_t splim2 = sqr(mhd->par.speedlimit / mhd->vvnorm);
  mrc_fld_data_t gamm   = mhd->par.gamm;
  mrc_fld_data_t d_i    = mhd->par.d_i;

  mrc_fld_data_t thx    = mhd->par.thx;
  mrc_fld_data_t eps    = 1e-9f;
  mrc_fld_data_t dt     = 1e10f;

  mrc_fld_data_t two_pi_d_i = 2. * M_PI * d_i;
  bool have_hall = d_i > 0.f;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);

    mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
      mrc_fld_data_t hh = 0.f;
      if (s_sw[0]) hh = mrc_fld_max(hh, PDE_INV_DX(ix));
      if (s_sw[1]) hh = mrc_fld_max(hh, PDE_INV_DY(iy));
      if (s_sw[2]) hh = mrc_fld_max(hh, PDE_INV_DZ(iz));

      mrc_fld_data_t rri = 1.f / RR_(x, ix,iy,iz, p);
      mrc_fld_data_t bb = (sqr(BTX_(x, ix,iy,iz, p)) + 
			   sqr(BTY_(x, ix,iy,iz, p)) +
			   sqr(BTZ_(x, ix,iy,iz, p)));
      mrc_fld_data_t rrvv = (sqr(RVX_(x, ix,iy,iz, p)) + 
			     sqr(RVY_(x, ix,iy,iz, p)) +
			     sqr(RVZ_(x, ix,iy,iz, p)));
      mrc_fld_data_t pp = (gamm - 1.f) * (EE_(x, ix,iy,iz, p) - .5f * rrvv * rri - .5f * bb);

      if (have_hall) {
	bb *= 1 + sqr(two_pi_d_i * hh);
      }      
      
      mrc_fld_data_t vv1 = mrc_fld_min(bb * rri, splim2);
      mrc_fld_data_t vv2 = gamm * pp * rri;
      mrc_fld_data_t vv = mrc_fld_sqrt(vv1 + vv2) + mrc_fld_sqrt(rrvv) * rri;

      dt = mrc_fld_min(dt, thx / mrc_fld_max(eps, hh*vv));
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  return dtn;
}

// ----------------------------------------------------------------------
// pde_mhd_get_dt_fcons_ct

// FIXME, take into account resistivity
// FIXME, rework

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_fcons_ct(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *Bcc)
{
  static const mrc_fld_data_t eps = 1.e-10; // FIXME

  mrc_fld_data_t qsq,asq,cf1sq,cf2sq,cf3sq;
  mrc_fld_data_t b1,b2,b3,bsq,tsum,tdif;
  mrc_fld_data_t max_v1=0.0,max_v2=0.0,max_v3=0.0,max_dti = 0.0;

  mrc_fld_data_t gamma = mhd->par.gamm;
  mrc_fld_data_t gamma_minus_1 = gamma - 1.;
  
  mrc_fld_data_t d_i = mhd->par.d_i;

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
					      - .5f*bsq), eps);
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

  // FOR diffusive portion of everything
  if (d_i > 0.0) {
    mrc_fld_data_t local_diff_dt, global_diff_dt;
    mrc_fld_data_t max_dti_diff = 0.0;

    for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
      double dx[3]; mrc_crds_get_dx(crds, p, dx);

      mrc_fld_foreach(x, i,j,k, 0, 0) {
	mrc_fld_data_t dxmin, dti_hall;

	dxmin = mrc_fld_min(dx[0], mrc_fld_min(dx[1], dx[2]));
	
	// FIXME: this could be streamlined by using the loop above, but
	//        keeping them separate is exactly analogous to athena
	/* Use maximum of face-centered fields (always larger than cell-centered B) */
	b1 = M3(Bcc, 0, i,j,k, p) + fabs(BX_(x, i,j,k, p) - M3(Bcc, 0, i,j,k, p));
	b2 = M3(Bcc, 1, i,j,k, p) + fabs(BY_(x, i,j,k, p) - M3(Bcc, 1, i,j,k, p));
	b3 = M3(Bcc, 2, i,j,k, p) + fabs(BZ_(x, i,j,k, p) - M3(Bcc, 2, i,j,k, p));
	bsq = sqr(b1) + sqr(b2) + sqr(b3);
	
	// FIXME? in athena it was (dxmin**2 / 6), but i don't think that's right?
	dti_hall = (d_i * bsq / RR_(x, i,j,k, p)) / (sqr(dxmin) / 16.0);
	max_dti_diff = mrc_fld_max(max_dti_diff, dti_hall);
      } mrc_fld_foreach_end;
    }
    local_diff_dt = cfl / max_dti_diff;
    MPI_Allreduce(&local_diff_dt, &global_diff_dt, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));
    global_dt = mrc_fld_min(global_dt, global_diff_dt);
  }

  return global_dt;
}

