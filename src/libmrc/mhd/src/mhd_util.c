
#include "mhd_util.h"

#include "ggcm_mhd_private.h"

#include <mrc_domain.h>
#include <mrc_crds.h>
#include <mrc_bits.h>
#include <mrc_fld_as_double.h> // FIXME

void
setup_mrc_fld_3d(struct mrc_fld *f, struct mrc_fld *f_tmpl, int nr_comps)
{
  mrc_fld_set_type(f, mrc_fld_type(f_tmpl));
  mrc_fld_set_param_obj(f, "domain", f_tmpl->_domain);
  mrc_fld_set_param_int(f, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(f, "nr_comps", nr_comps);
  mrc_fld_set_param_int(f, "nr_ghosts", f_tmpl->_nr_ghosts);
}

void
setup_mrc_fld_1d(struct mrc_fld *f, struct mrc_fld *f_tmpl, int nr_comps)
{
  int size = 0;
  for (int d = 0; d < 3; d++) {
    size = MAX(size, mrc_fld_spatial_dims(f_tmpl)[d] + 2 * f_tmpl->_nr_ghosts);
  }

  mrc_fld_set_type(f, mrc_fld_type(f_tmpl));
  mrc_fld_set_param_int_array(f , "dims", 2, (int []) { nr_comps, size });
  mrc_fld_set_param_int_array(f , "offs", 2, (int []) { 0, -f_tmpl->_nr_ghosts });
}

