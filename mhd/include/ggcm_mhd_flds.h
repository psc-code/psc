
#ifndef GGCM_MHD_FLDS_H
#define GGCM_MHD_FLDS_H

#include <mrc_obj.h>

// ======================================================================
// ggcm_mhd_flds
//
// This object stores the MHD fields (possibly wrapping the Fortran common block)

MRC_CLASS_DECLARE(ggcm_mhd_flds, struct ggcm_mhd_flds);

struct ggcm_mhd_flds *ggcm_mhd_flds_duplicate(struct ggcm_mhd_flds *flds);
void ggcm_mhd_flds_copy(struct ggcm_mhd_flds *to, struct ggcm_mhd_flds *from);
struct ggcm_mhd_flds *ggcm_mhd_flds_get_as(struct ggcm_mhd_flds *flds_base,
					   const char *type);
void ggcm_mhd_flds_put_as(struct ggcm_mhd_flds *flds,
			  struct ggcm_mhd_flds *flds_base);
struct mrc_f3 *ggcm_mhd_flds_get_mrc_f3(struct ggcm_mhd_flds *flds);
struct mrc_fld *ggcm_mhd_flds_get_mrc_fld(struct ggcm_mhd_flds *flds);

#endif

