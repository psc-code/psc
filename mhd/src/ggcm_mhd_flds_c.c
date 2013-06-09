
#include "ggcm_mhd_flds_private.h"

#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_io.h>
#include <assert.h>
#include <string.h>

// FIXME, duplicated
static const char *fldname[_NR_FLDS] = {
  [_RR1 ] = "rr1",
  [_RV1X] = "rv1x",
  [_RV1Y] = "rv1y",
  [_RV1Z] = "rv1z",
  [_UU1 ] = "uu1",
  [_B1X ] = "b1x",
  [_B1Y ] = "b1y",
  [_B1Z ] = "b1z",

  [_RR2 ] = "rr2",
  [_RV2X] = "rv2x",
  [_RV2Y] = "rv2y",
  [_RV2Z] = "rv2z",
  [_UU2 ] = "uu2",
  [_B2X ] = "b2x",
  [_B2Y ] = "b2y",
  [_B2Z ] = "b2z",

  [_YMASK] = "ymask",
  [_ZMASK] = "zmask",
  [_CMSV ] = "cmsv",

  [_RR  ] = "rr",
  [_PP  ] = "pp",
  [_VX  ] = "vx",
  [_VY  ] = "vy",
  [_VZ  ] = "vz",
  [_BX  ] = "bx",
  [_BY  ] = "by",
  [_BZ  ] = "bz",

  [_TMP1] = "tmp1",
  [_TMP2] = "tmp2",
  [_TMP3] = "tmp3",
  [_TMP4] = "tmp4",

  [_FLX ] = "ex",
  [_FLY ] = "ey",
  [_FLZ ] = "ez",

  [_CX  ] = "cx",
  [_CY  ] = "cy",
  [_CZ  ] = "cz",

  [_XTRA1] = "xtra1",
  [_XTRA2] = "xtra2",

  [_RESIS] = "resis",

  [_CURRX] = "currx",
  [_CURRY] = "curry",
  [_CURRZ] = "currz",

  [_RMASK] = "rmask",

  [_BDIPX] = "bdipx",
  [_BDIPY] = "bdipy",
  [_BDIPZ] = "bdipz",
};

// ======================================================================
// ggcm_mhd_flds subclass "c"

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_create

static void
ggcm_mhd_flds_c_create(struct ggcm_mhd_flds *flds)
{
  flds->fld = mrc_fld_create(ggcm_mhd_flds_comm(flds));
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_setup

static void
ggcm_mhd_flds_c_setup(struct ggcm_mhd_flds *flds)
{
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(flds->mhd->domain, &nr_patches);
  assert(nr_patches == 1);
  int *ldims = patches[0].ldims;
  mrc_fld_set_param_int_array(flds->fld, "dims", 4,
			      (int[4]) { ldims[0], ldims[1], ldims[2], _NR_FLDS });
  mrc_fld_set_param_int_array(flds->fld, "sw", 4,
			      (int[4]) { BND, BND, BND, 0 });
  flds->fld->_domain = flds->mhd->domain;
  mrc_fld_setup(flds->fld);
  for (int m = 0; m < _NR_FLDS; m++) {
    mrc_fld_set_comp_name(flds->fld, m, fldname[m]);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_destroy

static void
ggcm_mhd_flds_c_destroy(struct ggcm_mhd_flds *flds)
{
  mrc_fld_destroy(flds->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_write

static void
ggcm_mhd_flds_c_write(struct ggcm_mhd_flds *flds, struct mrc_io *io)
{
  mrc_io_write_ref(io, flds, "fld", flds->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_read

static void
ggcm_mhd_flds_c_read(struct ggcm_mhd_flds *flds, struct mrc_io *io)
{
  // FIXME, this calls out ::setup(), which we don't really want...
  ggcm_mhd_flds_read_super(flds, io);
  flds->fld = mrc_io_read_ref(io, flds, "fld", mrc_fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_copy

static void
ggcm_mhd_flds_c_copy(struct ggcm_mhd_flds *to, struct ggcm_mhd_flds *from_base)
{
  struct ggcm_mhd_flds *from = ggcm_mhd_flds_get_as(from_base, "c");
  mrc_fld_copy(to->fld, from->fld);
  ggcm_mhd_flds_put_as(from, from_base);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_copy_from_fortran

static void
ggcm_mhd_flds_c_copy_from_fortran(struct ggcm_mhd_flds *flds_c,
				  struct ggcm_mhd_flds *flds_fortran)
{
  mrc_fld_copy(flds_c->fld, ggcm_mhd_flds_get_mrc_fld(flds_fortran));
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_copy_to_fortran

static void
ggcm_mhd_flds_c_copy_to_fortran(struct ggcm_mhd_flds *flds_c,
				struct ggcm_mhd_flds *flds_fortran)
{
  mrc_fld_copy(ggcm_mhd_flds_get_mrc_fld(flds_fortran), flds_c->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds subclass "c"

static struct mrc_obj_method ggcm_mhd_flds_c_methods[] = {
  MRC_OBJ_METHOD("copy_to_fortran",   ggcm_mhd_flds_c_copy_to_fortran),
  MRC_OBJ_METHOD("copy_from_fortran", ggcm_mhd_flds_c_copy_from_fortran),
  {}
};

struct ggcm_mhd_flds_ops ggcm_mhd_flds_ops_c = {
  .name             = "c",
  .methods          = ggcm_mhd_flds_c_methods,
  .create           = ggcm_mhd_flds_c_create,
  .setup            = ggcm_mhd_flds_c_setup,
  .destroy          = ggcm_mhd_flds_c_destroy,
  .write            = ggcm_mhd_flds_c_write,
  .read             = ggcm_mhd_flds_c_read,
  .copy             = ggcm_mhd_flds_c_copy,
};
