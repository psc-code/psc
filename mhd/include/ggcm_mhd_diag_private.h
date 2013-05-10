
#ifndef GGCM_MHD_DIAG_PRIVATE_H
#define GGCM_MHD_DIAG_PRIVATE_H

#include "ggcm_mhd_diag.h"

#include <mrc_fld.h>

struct ggcm_mhd_diag {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_diag_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_diag);
  void (*run)(struct ggcm_mhd_diag *);
  void (*shutdown)(struct ggcm_mhd_diag *);
};

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;
extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_s2_ops;
extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_f2_ops;

void ggcm_mhd_diag_c_write_one_field(struct mrc_io *io, struct mrc_f3 *f, int m,
				     const char *name, float scale, int outtype,
				     float plane);
void ggcm_mhd_diag_c_write_one_f3(struct mrc_io *io, struct mrc_f3 *f,
				  int outtype, float plane);

#endif
