
// ======================================================================
// MHD parameters, we keep these around statically

static double s_gamma;  // adiabatic exponent

// ----------------------------------------------------------------------
// pde_mhd_setup

static void
pde_mhd_setup(struct ggcm_mhd *mhd)
{
  s_gamma = mhd->par.gamm;
}

