
#include "ggcm_mhd_flds_private.h"
#include "ggcm_mhd.h"

#include <mrc_io.h>
#include <mrc_profile.h>
#include <assert.h>
#include <string.h>

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
// ggcm_mhd_flds class

// ----------------------------------------------------------------------
// ggcm_mhd_flds_create

static void
_ggcm_mhd_flds_create(struct ggcm_mhd_flds *flds)
{
  flds->fld = mrc_fld_create(ggcm_mhd_flds_comm(flds));
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_setup

static void
_ggcm_mhd_flds_setup(struct ggcm_mhd_flds *flds)
{
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(flds->fld->_domain, &nr_patches);
  assert(nr_patches == 1);
  int *dims = patches[0].ldims;
  mrc_fld_set_param_int_array(flds->fld, "dims", 4,
			      (int[4]) { dims[0], dims[1], dims[2], _NR_FLDS });
  mrc_fld_set_param_int_array(flds->fld, "sw", 4,
			      (int[4]) { BND, BND, BND, 0 });

  for (int m = 0; m < _NR_FLDS; m++) {
    mrc_fld_set_comp_name(flds->fld, m, fldname[m]);
  }
  mrc_fld_setup(flds->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_destroy

static void
_ggcm_mhd_flds_destroy(struct ggcm_mhd_flds *flds)
{
  mrc_fld_destroy(flds->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_read

static void
_ggcm_mhd_flds_read(struct ggcm_mhd_flds *flds, struct mrc_io *io)
{
  flds->fld = mrc_io_read_ref(io, flds, "fld", mrc_fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_write

static void
_ggcm_mhd_flds_write(struct ggcm_mhd_flds *flds, struct mrc_io *io)
{
  mrc_io_write_ref(io, flds, "fld", flds->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_get_mrc_fld
//
// returns the underlying mrc_fld that contains the MHD field data
// only actually works (and should be used) if we know that these flds
// are of type "fortran" or "c"

struct mrc_fld *
ggcm_mhd_flds_get_mrc_fld(struct ggcm_mhd_flds *flds)
{
  assert(flds->fld);
  return flds->fld;
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_get_as
//
// convert flds_base to ggcm_mhd_flds of type "type"

struct ggcm_mhd_flds *
ggcm_mhd_flds_get_as(struct ggcm_mhd_flds *flds_base, const char *type)
{
  const char *type_base = ggcm_mhd_flds_type(flds_base);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return flds_base;

  // special case: when getting a "c" field from a "fortran" base,
  // we just copy the array pointer
  if (strcmp(type_base, "fortran") == 0 && strcmp(type, "c") == 0) {
    struct ggcm_mhd_flds *flds = ggcm_mhd_flds_create(ggcm_mhd_flds_comm(flds_base));
    ggcm_mhd_flds_set_type(flds, type);
    flds->fld->_domain = flds_base->fld->_domain;
    mrc_fld_set_array(flds->fld, flds_base->fld->_arr);
    ggcm_mhd_flds_setup(flds);
    return flds;
  }

  static int pr;
  if (!pr) {
    pr = prof_register("ggcm_mhd_flds_get_as", 1., 0, 0);
  }
  prof_start(pr);

  struct ggcm_mhd_flds *flds = ggcm_mhd_flds_create(ggcm_mhd_flds_comm(flds_base));
  ggcm_mhd_flds_set_type(flds, type);
  flds->fld->_domain = flds_base->fld->_domain;
  ggcm_mhd_flds_setup(flds);

  char s[strlen(type) + 12]; sprintf(s, "copy_to_%s", type);
  mrc_fld_copy_to_func_t copy_to = (mrc_fld_copy_to_func_t)
    mrc_fld_get_method(flds_base->fld, s);
  if (copy_to) {
    copy_to(flds_base->fld, flds->fld);
  } else {
    sprintf(s, "copy_from_%s", type_base);
    mrc_fld_copy_from_func_t copy_from = (mrc_fld_copy_from_func_t)
      mrc_fld_get_method(flds->fld, s);
    if (copy_from) {
      copy_from(flds->fld, flds_base->fld);
    } else {
      fprintf(stderr, "ERROR: no 'copy_to_%s' in ggcm_mhd_flds '%s' and "
	      "no 'copy_from_%s' in '%s'!\n",
	      type, ggcm_mhd_flds_type(flds_base), type_base, ggcm_mhd_flds_type(flds));
      assert(0);
    }
  }

  prof_stop(pr);
  return flds;
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_put_as
//
// after being done with the fields gotten from get_as(), need to put them
// back using this routine, which will copy the contents from flds back
// to flds_base

void
ggcm_mhd_flds_put_as(struct ggcm_mhd_flds *flds, struct ggcm_mhd_flds *flds_base)
{
  const char *type_base = ggcm_mhd_flds_type(flds_base);
  const char *type = ggcm_mhd_flds_type(flds);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return;

  // special case: when we originall got a "c" field from a "fortran" base,
  // we just copied the array pointer, so there's nothing to copy back
  if (strcmp(type_base, "fortran") == 0 && strcmp(type, "c") == 0) {
    ggcm_mhd_flds_destroy(flds);
    return;
  }

  static int pr;
  if (!pr) {
    pr = prof_register("ggcm_mhd_flds_put_as", 1., 0, 0);
  }
  prof_start(pr);

  char s[strlen(type) + 12]; sprintf(s, "copy_from_%s", type);
  mrc_fld_copy_from_func_t copy_from = (mrc_fld_copy_from_func_t)
    mrc_fld_get_method(flds_base->fld, s);
  if (copy_from) {
    copy_from(flds_base->fld, flds->fld);
  } else {
    sprintf(s, "copy_to_%s", type_base);
    mrc_fld_copy_to_func_t copy_to = (mrc_fld_copy_to_func_t)
      mrc_fld_get_method(flds->fld, s);
    if (copy_to) {
      copy_to(flds->fld, flds_base->fld);
    } else {
      fprintf(stderr, "ERROR: no 'copy_from_%s' in ggcm_mhd_flds '%s' and "
	      "no 'copy_to_%s' in '%s'!\n",
	      type, ggcm_mhd_flds_type(flds_base), type_base, ggcm_mhd_flds_type(flds));
      assert(0);
    }
  }

  ggcm_mhd_flds_destroy(flds);

  prof_stop(pr);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds_init

static void
ggcm_mhd_flds_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_flds, &ggcm_mhd_flds_ops_c);
}

// ----------------------------------------------------------------------
// ggcm_mhd_flds description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_flds, x)
static struct param ggcm_mhd_flds_descr[] = {
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_flds class description

struct mrc_class_ggcm_mhd_flds mrc_class_ggcm_mhd_flds = {
  .name             = "ggcm_mhd_flds",
  .size             = sizeof(struct ggcm_mhd_flds),
  .param_descr      = ggcm_mhd_flds_descr,
  .init             = ggcm_mhd_flds_init,
  .create           = _ggcm_mhd_flds_create,
  .setup            = _ggcm_mhd_flds_setup,
  .destroy          = _ggcm_mhd_flds_destroy,
  .read             = _ggcm_mhd_flds_read,
  .write            = _ggcm_mhd_flds_write,
};

