
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_vec.h>
#include <mrc_io.h>
#include <assert.h>
#include <string.h>

// ======================================================================
// mrc_fld subclass "mhd_pr_float"
//
// This subclass maintains mhd fields in primitive (pr) form
// Note that B is still staggered as usual, and the field indices start at
// 0, not the higher index that openggcm uses

struct mrc_fld_mhd_pr_float {
  struct ggcm_mhd *mhd;
};

#define mrc_fld_mhd_pr_float(fld) mrc_to_subobj(fld, struct mrc_fld_mhd_pr_float)

// ----------------------------------------------------------------------
// mrc_fld_mhd_pr_float_create

static void
mrc_fld_mhd_pr_float_create(struct mrc_fld *fld)
{
  fld->_data_type = MRC_NT_FLOAT;
  fld->_size_of_type = sizeof(float);
}

// ----------------------------------------------------------------------
// mrc_fld_mhd_pr_float_copy_from_float

static void
mrc_fld_mhd_pr_float_copy_from_float(struct mrc_fld *fld_pr, struct mrc_fld *fld_sc)
{
  return;
  struct mrc_fld *pr = mrc_fld_get_as(fld_pr, "mhd_pr_float");
  struct mrc_fld *sc = mrc_fld_get_as(fld_sc, "float");

  struct ggcm_mhd *mhd = mrc_fld_get_param_obj(fld_sc, "mhd");
  assert(mhd);
  float gamma = mhd->par.gamm;

  // FIXME, actually does "fc" for now
  mrc_fld_foreach(pr, ix, iy, iz, 1, 1) {
    RR1(pr, ix,iy,iz) = RR1 (sc, ix,iy,iz);
    V1X(pr, ix,iy,iz) = RV1X(sc, ix,iy,iz) / RR1(sc, ix,iy,iz);
    V1Y(pr, ix,iy,iz) = RV1Y(sc, ix,iy,iz) / RR1(sc, ix,iy,iz);
    V1Z(pr, ix,iy,iz) = RV1Z(sc, ix,iy,iz) / RR1(sc, ix,iy,iz);
    PP1(pr, ix,iy,iz) = (gamma - 1.f) * 
      (UU1 (sc, ix,iy,iz) -
      .5f * RR1(pr, ix, iy, iz) * (sqr(V1X(pr, ix,iy,iz)) +
				   sqr(V1Y(pr, ix,iy,iz)) +
				   sqr(V1Z(pr, ix,iy,iz))) -
      .5f * (sqr(.5*(B1X(sc, ix,iy,iz) + B1X(sc, ix+1,iy,iz))) +
	     sqr(.5*(B1Y(sc, ix,iy,iz) + B1Y(sc, ix,iy+1,iz))) +
	     sqr(.5*(B1Z(sc, ix,iy,iz) + B1Z(sc, ix,iy,iz+1)))));
    B1X (pr, ix,iy,iz) = B1X (sc, ix,iy,iz);
    B1Y (pr, ix,iy,iz) = B1Y (sc, ix,iy,iz);
    B1Z (pr, ix,iy,iz) = B1Z (sc, ix,iy,iz);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(pr, fld_pr);
  mrc_fld_put_as(sc, fld_sc);
}

// ----------------------------------------------------------------------
// mrc_fld_mhd_pr_float_copy_to_float

static void
mrc_fld_mhd_pr_float_copy_to_float(struct mrc_fld *fld_pr, struct mrc_fld *fld_sc)
{
  struct mrc_fld *pr = mrc_fld_get_as(fld_pr, "mhd_pr_float");
  struct mrc_fld *sc = mrc_fld_get_as(fld_sc, "float");

  struct ggcm_mhd *mhd = mrc_fld_get_param_obj(fld_sc, "mhd");
  assert(mhd);
  float gamma = mhd->par.gamm;

  // FIXME, actually does "fc" for now
  mrc_fld_foreach(sc, ix, iy, iz, 2, 2) {
    RR1 (sc, ix,iy,iz) = RR1(pr, ix,iy,iz);
    RV1X(sc, ix,iy,iz) = RR1(pr, ix,iy,iz) * V1X(pr, ix,iy,iz);
    RV1Y(sc, ix,iy,iz) = RR1(pr, ix,iy,iz) * V1Y(pr, ix,iy,iz);
    RV1Z(sc, ix,iy,iz) = RR1(pr, ix,iy,iz) * V1Z(pr, ix,iy,iz);
    UU1 (sc, ix,iy,iz) = PP1(pr, ix,iy,iz) / (gamma - 1.f) + 	
      .5f * RR1(pr, ix, iy, iz) * (sqr(V1X(pr, ix,iy,iz)) +
				   sqr(V1Y(pr, ix,iy,iz)) +
				   sqr(V1Z(pr, ix,iy,iz))) +
      .5f * (sqr(.5*(B1X(pr, ix,iy,iz) + B1X(pr, ix+1,iy,iz))) +
	     sqr(.5*(B1Y(pr, ix,iy,iz) + B1Y(pr, ix,iy+1,iz))) +
	     sqr(.5*(B1Z(pr, ix,iy,iz) + B1Z(pr, ix,iy,iz+1))));

    B1X (sc, ix,iy,iz) = B1X(pr, ix,iy,iz);
    B1Y (sc, ix,iy,iz) = B1Y(pr, ix,iy,iz);
    B1Z (sc, ix,iy,iz) = B1Z(pr, ix,iy,iz);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(pr, fld_pr);
  mrc_fld_put_as(sc, fld_sc);
}

// ----------------------------------------------------------------------
// mrc_fld subclass "mhd_pr_float" 

static struct mrc_obj_method mrc_fld_mhd_pr_float_methods[] = {
  MRC_OBJ_METHOD("copy_to_float",   mrc_fld_mhd_pr_float_copy_to_float),
  MRC_OBJ_METHOD("copy_from_float", mrc_fld_mhd_pr_float_copy_from_float),
  {}
};

#define VAR(x) (void *)offsetof(struct mrc_fld_mhd_pr_float, x)
static struct param mrc_fld_mhd_pr_float_descr[] = {
  { "mhd"             , VAR(mhd)          , PARAM_OBJ(ggcm_mhd)      },
  {},
};
#undef VAR

struct mrc_fld_ops mrc_fld_ops_mhd_pr_float = {
  .name             = "mhd_pr_float",
  .size             = sizeof(struct mrc_fld_mhd_pr_float),
  .param_descr      = mrc_fld_mhd_pr_float_descr,
  .methods          = mrc_fld_mhd_pr_float_methods,
  .create           = mrc_fld_mhd_pr_float_create,
  .vec_type         = "float",
};

