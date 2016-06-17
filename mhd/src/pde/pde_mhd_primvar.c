
// ----------------------------------------------------------------------
// patch_primvar

static void
patch_primvar(fld3d_t p_f, int m)
{
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;

  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_f,_RR, i,j,k) = F3S(p_f, m + _RR1, i,j,k);
    mrc_fld_data_t rri = 1.f / F3S(p_f, m + _RR1, i,j,k);
    F3S(p_f,_VX, i,j,k) = rri * F3S(p_f, m + _RV1X, i,j,k);
    F3S(p_f,_VY, i,j,k) = rri * F3S(p_f, m + _RV1Y, i,j,k);
    F3S(p_f,_VZ, i,j,k) = rri * F3S(p_f, m + _RV1Z, i,j,k);
    mrc_fld_data_t rvv =
      F3S(p_f,_VX, i,j,k) * F3S(p_f, m + _RV1X, i,j,k) +
      F3S(p_f,_VY, i,j,k) * F3S(p_f, m + _RV1Y, i,j,k) +
      F3S(p_f,_VZ, i,j,k) * F3S(p_f, m + _RV1Z, i,j,k);
    F3S(p_f,_PP, i,j,k) = gamma_m1 * (F3S(p_f, m + _UU1, i,j,k) - .5f * rvv);
    mrc_fld_data_t cs2 = mrc_fld_max(s_gamma * F3S(p_f,_PP, i,j,k) * rri, 0.f);
    F3S(p_f,_CMSV, i,j,k) = mrc_fld_sqrt(rvv * rri) + mrc_fld_sqrt(cs2);
  } fld3d_foreach_end;
}

