
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

  // dtn is the global new timestep
  if (dtn <= mhd->par.dtmin) {
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin. Dying now!\n");
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt %g -> %g, dtmin = %g\n",
	       mhd->dt, dtn, mhd->par.dtmin);
    ggcm_mhd_wrongful_death(mhd, mhd->fld, -1);
  }

  return dtn;
}

// ----------------------------------------------------------------------
// pde_mhd_get_dt_fcons

static mrc_fld_data_t _mrc_unused
pde_mhd_get_dt_fcons(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *zmask, 
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

  // dtn is the global new timestep
  if (dtn <= mhd->par.dtmin) {
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin. Dying now!\n");
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt %g -> %g, dtmin = %g\n",
	       mhd->dt, dtn, mhd->par.dtmin);
    ggcm_mhd_wrongful_death(mhd, mhd->fld, -1);
  }

  return dtn;
}

