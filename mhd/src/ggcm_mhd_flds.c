
#include "ggcm_mhd_flds_private.h"
#include "ggcm_mhd.h"

#include <mrc_io.h>
#include <mrc_profile.h>
#include <assert.h>
#include <string.h>

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

