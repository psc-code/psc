
#include "ggcm_mhd_flds_private.h"

#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_io.h>
#include <assert.h>
#include <string.h>

// ======================================================================
// ggcm_mhd_flds subclass "c"

// ----------------------------------------------------------------------
// ggcm_mhd_flds_c_create

static void
ggcm_mhd_flds_c_create(struct ggcm_mhd_flds *flds)
{
  mrc_fld_set_type(flds->fld, "float");
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
  .read             = ggcm_mhd_flds_c_read,
  .copy             = ggcm_mhd_flds_c_copy,
};
