
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"

#include <mrc_domain.h>
#include <mrc_fld_as_double.h>

#include <string.h>

// ----------------------------------------------------------------------
// ggcm_mhd_convert_sc_ggcm_from_primitive
//
// converts from primitive variables to semi-conservative in-place.
// No ghost points are set.

static void
ggcm_mhd_convert_sc_ggcm_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  mrc_fld_data_t gamma_m1 = mhd->par.gamm - 1.;

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  mrc_fld_foreach(fld, ix,iy,iz, 0, 1) {
    RVX(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VX(fld, ix,iy,iz);
    RVY(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VY(fld, ix,iy,iz);
    RVZ(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VZ(fld, ix,iy,iz);
    UU (fld, ix,iy,iz) = PP(fld, ix,iy,iz) / gamma_m1
      + .5*(sqr(RVX(fld, ix,iy,iz)) +
	    sqr(RVY(fld, ix,iy,iz)) +
	    sqr(RVZ(fld, ix,iy,iz))) / RR(fld, ix,iy,iz);
    BX(fld, ix-1,iy,iz) = BX(fld, ix,iy,iz);
    BY(fld, ix,iy-1,iz) = BY(fld, ix,iy,iz);
    BZ(fld, ix,iy,iz-1) = BZ(fld, ix,iy,iz);
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_sc_from_primitive
//
// converts from primitive variables to semi-conservative alt B in-place.
// No ghost points are set.

static void
ggcm_mhd_convert_sc_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  mrc_fld_data_t gamma_m1 = mhd->par.gamm - 1.;

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
    RVX(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VX(fld, ix,iy,iz);
    RVY(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VY(fld, ix,iy,iz);
    RVZ(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VZ(fld, ix,iy,iz);
    UU (fld, ix,iy,iz) = PP(fld, ix,iy,iz) / gamma_m1
      + .5*(sqr(RVX(fld, ix,iy,iz)) +
	    sqr(RVY(fld, ix,iy,iz)) +
	    sqr(RVZ(fld, ix,iy,iz))) / RR(fld, ix,iy,iz);
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_fc_from_primitive
//
// converts from primitive variables to fully-conservative in-place.
// No ghost points are set, the staggered B fields need to exist on all faces
// (that means one more than cell-centered dims)

static void
ggcm_mhd_convert_fc_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  mrc_fld_data_t gamma_m1 = mhd->par.gamm - 1.;

  struct mrc_fld *fld = mrc_fld_get_as(fld_base, FLD_TYPE);

  // don't go into ghost cells in invariant directions
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int dx = (gdims[0] == 1) ? 0 : 1;
  int dy = (gdims[1] == 1) ? 0 : 1;
  int dz = (gdims[2] == 1) ? 0 : 1;

  mrc_fld_foreach(fld, ix,iy,iz, 0, 0) {
    RVX(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VX(fld, ix,iy,iz);
    RVY(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VY(fld, ix,iy,iz);
    RVZ(fld, ix,iy,iz) = RR(fld, ix,iy,iz) * VZ(fld, ix,iy,iz);
    UU (fld, ix,iy,iz) = PP(fld, ix,iy,iz) / gamma_m1
      + .5*(sqr(.5*(BX(fld, ix,iy,iz) + BX(fld, ix+dx,iy,iz))) +
	    sqr(.5*(BY(fld, ix,iy,iz) + BY(fld, ix,iy+dy,iz))) +
	    sqr(.5*(BZ(fld, ix,iy,iz) + BZ(fld, ix,iy,iz+dz))))
      + .5*(sqr(RVX(fld, ix,iy,iz)) +
	    sqr(RVY(fld, ix,iy,iz)) +
	    sqr(RVZ(fld, ix,iy,iz))) / RR(fld, ix,iy,iz);
  } mrc_fld_foreach_end;
  
  mrc_fld_put_as(fld, fld_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_convert_from_primitive
//
// converts from primitive variables to the appropriate fully-conservative /
// semi-conservative state vector.

void
ggcm_mhd_convert_from_primitive(struct ggcm_mhd *mhd, struct mrc_fld *fld_base)
{
  int mhd_type;
  mrc_fld_get_param_int(fld_base, "mhd_type", &mhd_type);

  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM) {
    return ggcm_mhd_convert_sc_ggcm_from_primitive(mhd, fld_base);
  } else if (mhd_type == MT_FULLY_CONSERVATIVE) {
    return ggcm_mhd_convert_fc_from_primitive(mhd, fld_base);
  } else if (mhd_type == MT_SEMI_CONSERVATIVE) {
    return ggcm_mhd_convert_sc_from_primitive(mhd, fld_base);
  } else {
    assert(0);
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// ggcm_mhd_get_fld_as_fortran
//
// get mhd::fld as Fortran common block field

struct mrc_fld *
ggcm_mhd_get_fld_as_fortran(struct mrc_fld *mhd_fld)
{
  int mhd_type;
  mrc_fld_get_param_int(mhd_fld, "mhd_type", &mhd_type);

  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM) {
    int nr_ghosts;
    mrc_fld_get_param_int(mhd_fld, "nr_ghosts", &nr_ghosts);
    assert(nr_ghosts == 2);
    if (strcmp(mrc_fld_type(mhd_fld), "float") == 0) {
      return mhd_fld;
    } else {
      assert(0);
    }
  }

  struct mrc_fld *fld = mrc_fld_create(mrc_fld_comm(mhd_fld));
  mrc_fld_set_type(fld, "float");
  mrc_fld_set_param_obj(fld, "domain", mhd_fld->_domain);
  mrc_fld_set_param_int(fld, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(fld, "nr_comps", _NR_FLDS);
  mrc_fld_set_param_int(fld, "nr_ghosts", 2);
  mrc_fld_setup(fld);

  if (mhd_type == MT_SEMI_CONSERVATIVE) {
    struct mrc_fld *f = mrc_fld_get_as(mhd_fld, "float");
    struct mrc_fld *ff = mrc_fld_get_as(fld, "float");
    for (int m = 0; m < _NR_FLDS; m++) {
      if (m == _B1X || m == _B2X) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(ff, m, ix,iy,iz) = MRC_F3(f, m, ix+1,iy,iz);
	} mrc_fld_foreach_end;
      } else if (m == _B1Y || m == _B2Y) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(ff, m, ix,iy,iz) = MRC_F3(f, m, ix,iy+1,iz);
	} mrc_fld_foreach_end;
      } else if (m == _B1Z || m == _B2Z) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(ff, m, ix,iy,iz) = MRC_F3(f, m, ix,iy,iz+1);
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(ff, m, ix,iy,iz) = MRC_F3(f, m, ix,iy,iz);
	} mrc_fld_foreach_end;
      }
    }
    mrc_fld_put_as(f, mhd_fld);
    mrc_fld_put_as(ff, fld);
  } else {
    assert(0);
  }

  return fld;
}

// ----------------------------------------------------------------------
// ggcm_mhd_put_fld_as_fortran

void
ggcm_mhd_put_fld_as_fortran(struct mrc_fld *mhd_fld, struct mrc_fld *fld)
{
  int mhd_type;
  mrc_fld_get_param_int(mhd_fld, "mhd_type", &mhd_type);

  if (mhd_type == MT_SEMI_CONSERVATIVE_GGCM) {
    if (strcmp(mrc_fld_type(mhd_fld), "float") == 0) {
      assert(fld == mhd_fld);
      return;
    } else {
      assert(0);
    }
  }

  if (mhd_type == MT_SEMI_CONSERVATIVE) {
    struct mrc_fld *f = mrc_fld_get_as(mhd_fld, "float");
    struct mrc_fld *ff = mrc_fld_get_as(fld, "float");
    for (int m = 0; m < _NR_FLDS; m++) {
      if (m == _B1X || m == _B2X) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(f, m, ix,iy,iz) = MRC_F3(ff, m, ix-1,iy,iz);
	} mrc_fld_foreach_end;
      } else if (m == _B1Y || m == _B2Y) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(f, m, ix,iy,iz) = MRC_F3(ff, m, ix,iy-1,iz);
	} mrc_fld_foreach_end;
      } else if (m == _B1Z || m == _B2Z) {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(f, m, ix,iy,iz) = MRC_F3(ff, m, ix,iy,iz-1);
	} mrc_fld_foreach_end;
      } else {
	mrc_fld_foreach(f, ix,iy,iz, 2, 2) {
	  MRC_F3(f, m, ix,iy,iz) = MRC_F3(ff, m, ix,iy,iz);
	} mrc_fld_foreach_end;
      }
    }
    mrc_fld_put_as(f, mhd_fld);
    mrc_fld_put_as(ff, fld);
  } else {
    assert(0);
  }

  mrc_fld_destroy(fld);
}

