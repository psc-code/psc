
#ifndef GGCM_MHD_CRDS_PRIVATE_H
#define GGCM_MHD_CRDS_PRIVATE_H

#include "ggcm_mhd_crds.h"

#include "ggcm_mhd_defs.h"

struct ggcm_mhd_crds {
  struct mrc_obj obj;

  bool legacy_fd1;
  struct mrc_domain *domain;

  struct mrc_fld *f1[3];
  struct mrc_fld *global_f1[3];

  // FIXME, this isn't generic for all crds types
  float *_crdx[NR_CRDS], *_crdy[NR_CRDS], *_crdz[NR_CRDS];
};

struct ggcm_mhd_crds_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_crds);
};

extern struct ggcm_mhd_crds_ops ggcm_mhd_crds_c_ops;

#endif
