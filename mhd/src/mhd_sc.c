
// ----------------------------------------------------------------------
// zmaskn

static void __unused
zmaskn(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  float va02i = 1.f / sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  float eps   = 1e-15f;

  mrc_fld_foreach(x, ix,iy,iz, 1, 1) {
    mrc_fld_data_t bb = (sqr(.5f * (BX(x, ix,iy,iz) + BX(x, ix-1,iy,iz))) +
			 sqr(.5f * (BY(x, ix,iy,iz) + BY(x, ix,iy-1,iz))) +
			 sqr(.5f * (BZ(x, ix,iy,iz) + BZ(x, ix,iy,iz-1))));
    float rrm = fmaxf(eps, bb * va02i);
    F3(x, _ZMASK, ix,iy,iz) = F3(x, _YMASK, ix,iy,iz) *
      fminf(1.f, RR(x, ix,iy,iz) / rrm);
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// newstep_sc

static mrc_fld_data_t __unused
newstep_sc(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  float *fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  float *fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  float *fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  mrc_fld_data_t splim2 = sqr(mhd->par.speedlimit / mhd->par.vvnorm);
  mrc_fld_data_t gamm   = mhd->par.gamm;
  mrc_fld_data_t thx    = mhd->par.thx;
  mrc_fld_data_t eps    = 1e-9f;
  mrc_fld_data_t dt     = 1e10f;

  mrc_fld_foreach(x, ix, iy, iz, 0, 0) {
    mrc_fld_data_t hh = mrc_fld_max(mrc_fld_max(fd1x[ix], fd1y[iy]), fd1z[iz]);
    mrc_fld_data_t rri = 1.f / mrc_fld_abs(RR(x, ix,iy,iz)); // FIXME abs necessary?
    mrc_fld_data_t bb = (sqr(.5f * (BX(x, ix,iy,iz) + BX(x, ix-1,iy,iz))) + 
			 sqr(.5f * (BY(x, ix,iy,iz) + BY(x, ix,iy-1,iz))) +
			 sqr(.5f * (BZ(x, ix,iy,iz) + BZ(x, ix,iy,iz-1))));
    mrc_fld_data_t vv1 = mrc_fld_min(bb * rri, splim2);
    mrc_fld_data_t rrvv = (sqr(RVX(x, ix,iy,iz)) + 
			   sqr(RVY(x, ix,iy,iz)) +
			   sqr(RVZ(x, ix,iy,iz)));
    mrc_fld_data_t pp = (gamm - 1.f) * (UU(x, ix,iy,iz) - .5f * rrvv * rri);
    mrc_fld_data_t vv2 = gamm * pp * rri;
    mrc_fld_data_t vv = mrc_fld_sqrt(vv1 + vv2) + mrc_fld_sqrt(rrvv) * rri;
    vv = mrc_fld_max(eps, vv);

    mrc_fld_data_t tt = thx / mrc_fld_max(eps, hh*vv*F3(x, _ZMASK, ix,iy,iz));
    dt = mrc_fld_min(dt, tt);
  } mrc_fld_foreach_end;

  mrc_fld_data_t dtn;
  MPI_Allreduce(&dt, &dtn, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));

  return dtn;
}

