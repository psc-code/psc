
// ----------------------------------------------------------------------
// patch_newstep_scons_ggcm

static mrc_fld_data_t
patch_newstep_scons_ggcm(fld3d_t p_f)
{
  mrc_fld_data_t dt = 1e10f;

  mrc_fld_data_t splim2   = sqr(s_speedlimit_code);
  mrc_fld_data_t isphere2 = sqr(s_isphere);
  mrc_fld_data_t va02i    = 1.f / splim2;
  mrc_fld_data_t eps      = 1e-9f;
  mrc_fld_data_t epsz     = 1e-15f;
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;

  mrc_fld_data_t two_pi_d_i = 2. * M_PI * s_d_i;
  bool have_hall = s_d_i > 0.f;

  fld3d_foreach(i, j, k, 0, 0) {
    mrc_fld_data_t hh = mrc_fld_max(mrc_fld_max(FD1X(i), FD1Y(j)), FD1Z(k));
    mrc_fld_data_t rri = 1.f / mrc_fld_abs(F3S(p_f, RR, i,j,k)); // FIME abs necessary?
    mrc_fld_data_t bb = 
      sqr(.5f * (F3S(p_f, BX, i,j,k) + F3S(p_f, BX, i-1,j,k))) + 
      sqr(.5f * (F3S(p_f, BY, i,j,k) + F3S(p_f, BY, i,j-1,k))) +
      sqr(.5f * (F3S(p_f, BZ, i,j,k) + F3S(p_f, BZ, i,j,k-1)));
    if (have_hall) {
      bb *= 1 + sqr(two_pi_d_i * hh);
    }
    mrc_fld_data_t vv1 = fminf(bb * rri, splim2);
    
    mrc_fld_data_t rv2 = 
      sqr(F3S(p_f, RVX, i,j,k)) + sqr(F3S(p_f, RVY, i,j,k)) + sqr(F3S(p_f, RVZ, i,j,k));
    mrc_fld_data_t rvv = rri * rv2;
    mrc_fld_data_t pp = gamma_m1 * (F3S(p_f, UU, i,j,k) - .5f * rvv);
    mrc_fld_data_t vv2 = s_gamma * mrc_fld_max(0.f, pp) * rri;
    mrc_fld_data_t vv3 = rri * sqrtf(rv2);
    mrc_fld_data_t vv = sqrtf(vv1 + vv2) + vv3;
    vv = mrc_fld_max(eps, vv);
    
    mrc_fld_data_t ymask = 1.f;
    if (FX2X(i) + FX2Y(j) + FX2Z(k) < isphere2)
      ymask = 0.f;
    
    mrc_fld_data_t rrm = mrc_fld_max(epsz, bb * va02i);
    mrc_fld_data_t zmask = ymask * fminf(1.f, F3S(p_f, RR, i,j,k) / rrm);
    
    mrc_fld_data_t tt = s_cfl / mrc_fld_max(eps, hh*vv*zmask);
    dt = mrc_fld_min(dt, tt);
  } fld3d_foreach_end;

  return dt;
}

