
#ifndef GGCM_MHD_FLDS_H
#define GGCM_MHD_FLDS_H

#include <mrc_obj.h>

// ======================================================================
// ggcm_mhd_flds
//
// This object stores the MHD fields (possibly wrapping the Fortran common block)

MRC_CLASS_DECLARE(ggcm_mhd_flds, struct ggcm_mhd_flds);

struct mrc_fld *ggcm_mhd_flds_get_mrc_fld(struct ggcm_mhd_flds *flds);

#endif

