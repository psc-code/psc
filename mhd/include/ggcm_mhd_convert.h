
#ifndef GGCM_MHD_CONVERT_H
#define GGCM_MHD_CONVERT_H

static mrc_fld_data_t cvt_gamma;
static mrc_fld_data_t cvt_gamma_m1;
static mrc_fld_data_t cvt_gamma_m1_inv;

static inline void ggcm_mhd_convert_setup(struct ggcm_mhd *mhd)
{
  cvt_gamma = mhd->par.gamm;
  cvt_gamma_m1 = cvt_gamma - 1.f;
  cvt_gamma_m1_inv = 1.f / cvt_gamma_m1;
}

static inline void
convert_state_from_prim_scons(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
  assert(cvt_gamma);
  
  state[RR ] = prim[RR];
  state[RVX] = prim[RR] * prim[VX];
  state[RVY] = prim[RR] * prim[VY];
  state[RVZ] = prim[RR] * prim[VZ];
  state[UU ] = prim[PP] * cvt_gamma_m1_inv +
    + .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]));
}
      
static inline void
convert_state_from_prim_fcons(mrc_fld_data_t state[8], mrc_fld_data_t prim[8])
{
  assert(cvt_gamma);

  state[RR ] = prim[RR];
  state[RVX] = prim[RR] * prim[VX];
  state[RVY] = prim[RR] * prim[VY];
  state[RVZ] = prim[RR] * prim[VZ];
  state[EE ] = prim[PP] * cvt_gamma_m1_inv
    + .5f * prim[RR] * (sqr(prim[VX]) + sqr(prim[VY]) + sqr(prim[VZ]))
    + .5f * (sqr(prim[BX]) + sqr(prim[BY]) + sqr(prim[BZ]));
}

#endif
