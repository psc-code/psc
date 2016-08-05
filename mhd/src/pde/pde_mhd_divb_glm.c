
#ifndef PDE_MHD_DIVB_GLM_C
#define PDE_MHD_DIVB_GLM_C

#include "pde/pde_mhd_convert.c"

// ======================================================================
// divb cleaning with generalized langrangian multiplier
// Dedner et al, JCP 2002

// ----------------------------------------------------------------------
// mhd_divb_glm_source
//
// solve source term in PSI equation
// Eq. (42) Dedner et al

static void _mrc_unused
mhd_divb_glm_source(struct mrc_fld *x, double dt)
{
  if (s_opt_divb != OPT_DIVB_GLM) {
    return;
  }

  //  mprintf("glm_ch %g\n", s_divb_glm_ch);

  mrc_fld_data_t dx = s_g_dxyzmin;
  // This is following Mignone et al, JCP 2010 instead of Dedner directly
  mrc_fld_data_t fac = exp(-dt / dx * s_divb_glm_ch * s_divb_glm_alpha);
  
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    mrc_fld_foreach(x, i,j,k, 0, 0) {
      M3(x, PSI, i,j,k, p) *= fac;
    } mrc_fld_foreach_end;
  }
}

// ----------------------------------------------------------------------
// mhd_divb_glm_riemann
//
// (exact) Riemann solve for the Bnormal / PSI part of the equations
// Eq. (42) Dedner et al

static void _mrc_unused
mhd_divb_glm_riemann(fld1d_state_t Ul, fld1d_state_t Ur, fld1d_state_t Wl, fld1d_state_t Wr, int ib, int ie)
{
  if (s_opt_divb != OPT_DIVB_GLM) {
    return;
  }

  for (int i = ib; i < ie; i++){
    mrc_fld_data_t Bm, psim;
    Bm =   (  .5f                 * (F1S(Wl, BX , i) + F1S(Wr, BX , i)) 
	    - .5f / s_divb_glm_ch * (F1S(Wr, PSI, i) - F1S(Wl, PSI, i)));
    psim = (  .5f                 * (F1S(Wl, PSI, i) + F1S(Wr, PSI, i))
	    - .5f * s_divb_glm_ch * (F1S(Wr, BX , i) - F1S(Wl, BX , i)));

    F1S(Wl, BX , i) = Bm;
    F1S(Wr, BX , i) = Bm;
    F1S(Wl, PSI, i) = psim;
    F1S(Wr, PSI, i) = psim;
  }

  mhd_cons_from_prim(Ul, Wl, ib, ie);
  mhd_cons_from_prim(Ur, Wr, ib, ie);
}

#endif
