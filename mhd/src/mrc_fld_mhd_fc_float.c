
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_vec.h>
#include <mrc_io.h>
#include <assert.h>
#include <string.h>

// ======================================================================
// mrc_fld subclass "mhd_fc_float"
//
// This subclass maintains mhd fields in fully-conservative (fc) form

struct mrc_fld_mhd_fc_float {
};

#define mrc_fld_mhd_fc_float(fld) mrc_to_subobj(fld, struct mrc_fld_mhd_fc_float)

// ----------------------------------------------------------------------
// mrc_fld_mhd_fc_float_create

static void
mrc_fld_mhd_fc_float_create(struct mrc_fld *fld)
{
  fld->_data_type = MRC_NT_FLOAT;
  fld->_size_of_type = sizeof(float);
}

// ----------------------------------------------------------------------
// mrc_fld_mhd_fc_float_copy_from_float

static void
mrc_fld_mhd_fc_float_copy_from_float(struct mrc_fld *fld_fc, struct mrc_fld *fld_sc)
{
  struct mrc_fld *fc = mrc_fld_get_as(fld_fc, "mhd_fc_float");
  struct mrc_fld *sc = mrc_fld_get_as(fld_sc, "float");

  mrc_fld_foreach(fc, ix, iy, iz, 2, 2) {
    RR1 (fc, ix,iy,iz) = RR1 (sc, ix,iy,iz);
    RV1X(fc, ix,iy,iz) = RV1X(sc, ix,iy,iz);
    RV1Y(fc, ix,iy,iz) = RV1Y(sc, ix,iy,iz);
    RV1Z(fc, ix,iy,iz) = RV1Z(sc, ix,iy,iz);
    UU1 (fc, ix,iy,iz) = UU1 (sc, ix,iy,iz) + 
      00*.5f * (sqr(.5*(B1X(sc, ix,iy,iz) + B1X(sc, ix+1,iy,iz))) +
	     sqr(.5*(B1Y(sc, ix,iy,iz) + B1Y(sc, ix,iy+1,iz))) +
	     sqr(.5*(B1Z(sc, ix,iy,iz) + B1Z(sc, ix,iy,iz+1))));
    B1X (fc, ix,iy,iz) = B1X (sc, ix,iy,iz);
    B1Y (fc, ix,iy,iz) = B1Y (sc, ix,iy,iz);
    B1Z (fc, ix,iy,iz) = B1Z (sc, ix,iy,iz);
    for (int m = 8; m < 11; m++) {
      MRC_F3(fc, m, ix,iy,iz) = MRC_F3(sc, m, ix,iy,iz);
    }
    MRC_F3(fc, _YMASK, ix,iy,iz) = MRC_F3(sc, _YMASK, ix,iy,iz);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fc, fld_fc);
  mrc_fld_put_as(sc, fld_sc);
}

// ----------------------------------------------------------------------
// mrc_fld_mhd_fc_float_copy_to_float

static void
mrc_fld_mhd_fc_float_copy_to_float(struct mrc_fld *fld_fc, struct mrc_fld *fld_sc)
{
  struct mrc_fld *fc = mrc_fld_get_as(fld_fc, "mhd_fc_float");
  struct mrc_fld *sc = mrc_fld_get_as(fld_sc, "float");

  mrc_fld_foreach(sc, ix, iy, iz, 2, 2) {
    RR1 (sc, ix,iy,iz) = RR1 (fc, ix,iy,iz);
    RV1X(sc, ix,iy,iz) = RV1X(fc, ix,iy,iz);
    RV1Y(sc, ix,iy,iz) = RV1Y(fc, ix,iy,iz);
    RV1Z(sc, ix,iy,iz) = RV1Z(fc, ix,iy,iz);
    UU1 (sc, ix,iy,iz) = UU1 (fc, ix,iy,iz) -
      00*.5f * (sqr(.5*(B1X(sc, ix,iy,iz) + B1X(fc, ix+1,iy,iz))) +
	     sqr(.5*(B1Y(sc, ix,iy,iz) + B1Y(fc, ix,iy+1,iz))) +
	     sqr(.5*(B1Z(sc, ix,iy,iz) + B1Z(fc, ix,iy,iz+1))));
    B1X (sc, ix,iy,iz) = B1X (fc, ix,iy,iz);
    B1Y (sc, ix,iy,iz) = B1Y (fc, ix,iy,iz);
    B1Z (sc, ix,iy,iz) = B1Z (fc, ix,iy,iz);
    for (int m = 8; m < 11; m++) {
      MRC_F3(sc, m, ix,iy,iz) = MRC_F3(fc, m, ix,iy,iz);
    }
    MRC_F3(sc, _YMASK, ix,iy,iz) = MRC_F3(fc, _YMASK, ix,iy,iz);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(fc, fld_fc);
  mrc_fld_put_as(sc, fld_sc);
}

// ----------------------------------------------------------------------
// mrc_fld subclass "mhd_fc_float" 

static struct mrc_obj_method mrc_fld_mhd_fc_float_methods[] = {
  MRC_OBJ_METHOD("copy_to_float",   mrc_fld_mhd_fc_float_copy_to_float),
  MRC_OBJ_METHOD("copy_from_float", mrc_fld_mhd_fc_float_copy_from_float),
  {}
};

struct mrc_fld_ops mrc_fld_ops_mhd_fc_float = {
  .name             = "mhd_fc_float",
  .size             = sizeof(struct mrc_fld_mhd_fc_float),
  .methods          = mrc_fld_mhd_fc_float_methods,
  .create           = mrc_fld_mhd_fc_float_create,
  .vec_type         = "float",
};

