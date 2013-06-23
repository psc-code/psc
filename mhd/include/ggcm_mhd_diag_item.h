
#ifndef GGCM_MHD_DIAG_ITEM_H
#define GGCM_MHD_DIAG_ITEM_H

#include <mrc_obj.h>
#include <mrc_fld.h>

// ======================================================================
// ggcm_mhd_diag_item

MRC_CLASS_DECLARE(ggcm_mhd_diag_item, struct ggcm_mhd_diag_item);

void ggcm_mhd_diag_item_run(struct ggcm_mhd_diag_item *item,
			    struct mrc_io *io, struct mrc_fld *f,
			    int diag_type, float plane);

#endif
