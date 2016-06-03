
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
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;
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
      mrc_fld_data_t pp = gamma_m1 * (EE_(x, ix,iy,iz, p) - .5f * rrvv * rri - .5f * bb);

      if (have_hall) {
	bb *= 1 + sqr(two_pi_d_i * hh);
      }      
      
      mrc_fld_data_t vv1 = mrc_fld_min(bb * rri, splim2);
      mrc_fld_data_t cs2 = s_gamma * pp * rri;
      mrc_fld_data_t vv = mrc_fld_sqrt(vv1 + cs2) + mrc_fld_sqrt(rrvv) * rri;

      dt = mrc_fld_min(dt, thx / mrc_fld_max(eps, hh*vv));
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

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
