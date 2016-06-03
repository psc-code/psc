
static inline mrc_fld_data_t
newstep_sc_inl(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *zmask, 
	       int m_zmask, struct mrc_fld *b0)
{
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
      mrc_fld_data_t bb = (sqr(.5f * (BT(x, 0, ix,iy,iz, p) + BT(x, 0, ix+dx,iy,iz, p))) + 
			   sqr(.5f * (BT(x, 1, ix,iy,iz, p) + BT(x, 1, ix,iy+dy,iz, p))) +
			   sqr(.5f * (BT(x, 2, ix,iy,iz, p) + BT(x, 2, ix,iy,iz+dz, p))));
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

  mrc_fld_data_t dtmin  = mhd->par.dtmin;
  // dtn is the global new timestep
  if (dtn <= dtmin) {
    // dt is local new timestep
    if (dt <= dtmin) {
      char filename[20];
      int rank;
      FILE *file;

      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      sprintf(filename, "dtmin.%d", rank);
      file = fopen(filename, "w");

      fprintf(file, "istep,dtmin %d %f\n",mhd->istep,dtmin);
      for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
	float *fd1x = ggcm_mhd_crds_get_crd_p(mhd->crds, 0, FD1, p);
	float *fd1y = ggcm_mhd_crds_get_crd_p(mhd->crds, 1, FD1, p);
	float *fd1z = ggcm_mhd_crds_get_crd_p(mhd->crds, 2, FD1, p);

	// FIXME: this re-implements all the math above
	mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
	  mrc_fld_data_t hh = fmaxf(fmaxf(fd1x[ix], fd1y[iy]), fd1z[iz]);
	  // FIXME, sqrtf() is not nec, we square the result again
	  mrc_fld_data_t vv =
	    mrc_fld_sqrt((sqr(BX(x, ix,iy,iz)) + sqr(BY(x, ix,iy,iz)) + sqr(BZ(x, ix,iy,iz)) +
			  gamm * fmaxf(0.f, PP(x, ix,iy,iz))) /
			 fabsf(RR(x, ix,iy,iz))) +
	    mrc_fld_sqrt(sqr(VX(x, ix,iy,iz)) + sqr(VY(x, ix,iy,iz)) + sqr(VZ(x, ix,iy,iz)));
	  vv = fmaxf(eps, vv);
	  mrc_fld_data_t zm = 1.f;
	  if (zmask) {
	    zm = F3(zmask, m_zmask, ix,iy,iz);
	  }
	  mrc_fld_data_t tt = thx / fmaxf(eps, hh*vv*zm);
	  if (tt <= dtmin) {
	    fprintf(file,"%s %7d %7d %7d\n"," dtmin at ", ix, iy, ix);
	    fprintf(file,"%s %12.5g %12.5g %12.5g\n"," hh,vv,tt ", hh, vv, tt);
	    fprintf(file,"%s %12.5g %12.5g %12.5g\n"," bx,y,z   ", BX(x, ix,iy,iz)*mhd->bbnorm, BY(x, ix,iy,iz)*mhd->bbnorm, BZ(x, ix,iy,iz)*mhd->bbnorm);
	    fprintf(file,"%s %12.5g %12.5g %12.5g\n"," vx,y,z   ", VX(x, ix,iy,iz)*mhd->vvnorm, VY(x, ix,iy,iz)*mhd->vvnorm, VZ(x, ix,iy,iz)*mhd->vvnorm);
	    fprintf(file,"%s %12.5g %12.5g %12.5g\n"," rr,pp,zm ", RR(x, ix,iy,iz)*mhd->rrnorm, PP(x, ix,iy,iz)*mhd->ppnorm, zm);
	    float *fx1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FX1);
	    float *fx1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FX1);
	    float *fx1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FX1);
	    fprintf(file,"%s %12.5g %12.5g %12.5g %12.5g\n"," xx,yy,zz,r ",fx1x[ix],fx1y[iy],fx1z[iz],sqrtf(sqr(fx1x[ix])+sqr(fx1y[iy])+sqr(fx1z[iz])));
	  }
	} mrc_fld_foreach_end;
      }
    }
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin. Dying now!\n");
    mpi_printf(ggcm_mhd_comm(mhd), "!!! dt %g -> %g, dtmin = %g\n",
	       mhd->dt, dtn, dtmin);
    ggcm_mhd_wrongful_death(mhd, mhd->fld, -1);
  }

  return dtn;
}

static mrc_fld_data_t _mrc_unused
newstep_sc(struct ggcm_mhd *mhd, struct mrc_fld *x, struct mrc_fld *zmask, 
	   int m_zmask)
{
  if (mhd->b0) {
    return newstep_sc_inl(mhd, x, zmask, m_zmask, mhd->b0);
  } else {
    return newstep_sc_inl(mhd, x, zmask, m_zmask, NULL);
  }
}

