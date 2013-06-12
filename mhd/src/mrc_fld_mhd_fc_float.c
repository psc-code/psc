
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
// mrc_fld_mhd_fc_float_setup

static void
mrc_fld_mhd_fc_float_setup(struct mrc_fld *fld)
{
  mrc_fld_setup_super(fld);

  mrc_vec_set_type(fld->_vec, "float");
  mrc_vec_set_param_int(fld->_vec, "len", fld->_len);
  mrc_fld_setup_member_objs(fld); // sets up our .vec member

  fld->_arr = mrc_vec_get_array(fld->_vec);
}

// ----------------------------------------------------------------------
// mrc_fld_mhd_fc_float_copy_from_float

static void
mrc_fld_mhd_fc_float_copy_from_float(struct mrc_fld *fld_mhd_fc_float,
				     struct mrc_fld *fld_float)
{
  mrc_fld_copy(fld_mhd_fc_float, fld_float);
}

// ----------------------------------------------------------------------
// mrc_fld_mhd_fc_float_copy_to_float

static void
mrc_fld_mhd_fc_float_copy_to_float(struct mrc_fld *fld_mhd_fc_float,
				   struct mrc_fld *fld_float)
{
  mrc_fld_copy(fld_float, fld_mhd_fc_float);
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
  .setup            = mrc_fld_mhd_fc_float_setup,
};

