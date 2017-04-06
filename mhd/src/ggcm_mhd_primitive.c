
#include "ggcm_mhd_private.h"

#include <mrc_fld_as_double.h>

// ----------------------------------------------------------------------
// ggcm_mhd_calc_rr

void
ggcm_mhd_calc_rr(struct ggcm_mhd *mhd, struct mrc_fld *rr_base,
		 struct mrc_fld *fld_base)
{
  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);

  struct mrc_fld *rr = mrc_fld_get_as(rr_base, FLD_TYPE);
  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    ggcm_mhd_calc_rr_scons(mhd, rr, fld);
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      ggcm_mhd_calc_rr_fcons_fc(mhd, rr, fld);
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      ggcm_mhd_calc_rr_fcons_cc(mhd, rr, fld);
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    ggcm_mhd_calc_rr_gkeyll(mhd, rr, fld);
  } else {
    assert(0);
  }

  mrc_fld_put_as(rr, rr_base);
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_v

void
ggcm_mhd_calc_v(struct ggcm_mhd *mhd, struct mrc_fld *v_base,
		struct mrc_fld *fld_base)
{
  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);

  struct mrc_fld *v = mrc_fld_get_as(v_base, FLD_TYPE);
  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    ggcm_mhd_calc_v_scons(mhd, v, fld);
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      ggcm_mhd_calc_v_fcons_fc(mhd, v, fld);
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      ggcm_mhd_calc_v_fcons_cc(mhd, v, fld);
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    ggcm_mhd_calc_v_gkeyll(mhd, v, fld);
  } else {
    assert(0);
  }

  mrc_fld_put_as(v, v_base);
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_calc_pp

void
ggcm_mhd_calc_pp(struct ggcm_mhd *mhd, struct mrc_fld *pp_base,
		 struct mrc_fld *fld_base)
{
  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);

  struct mrc_fld *pp = mrc_fld_get_as(pp_base, FLD_TYPE);
  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  if (MT_FORMULATION(mhd_type) == MT_FORMULATION_SCONS) {
    ggcm_mhd_calc_pp_scons(mhd, pp, fld);
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_FCONS) {
    if (MT_BGRID(mhd_type) == MT_BGRID_FC) {
      ggcm_mhd_calc_pp_fcons_fc(mhd, pp, fld);
    } else if (MT_BGRID(mhd_type) == MT_BGRID_CC) {
      ggcm_mhd_calc_pp_fcons_cc(mhd, pp, fld);
    }
  } else if (MT_FORMULATION(mhd_type) == MT_FORMULATION_GKEYLL) {
    ggcm_mhd_calc_pp_gkeyll(mhd, pp, fld);
  } else {
    assert(0);
  }

  mrc_fld_put_as(pp, pp_base);
  mrc_fld_put_as(fld, fld_base);
}

