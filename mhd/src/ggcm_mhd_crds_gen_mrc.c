
#include "ggcm_mhd_crds_gen_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds_private.h"

#include <mrc_domain.h>
#include <assert.h>

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_mrc_run
//
// initializes the Fortran common block coord arrays FX1, FD1 from
// using C grid generation

static void
ggcm_mhd_crds_gen_mrc_run(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(crds->mhd->domain, 0, &info);
  struct mrc_crds *mrc_crds = mrc_domain_get_crds(crds->mhd->domain);
  int *ldims = info.ldims;

  float *fxx1 = ggcm_mhd_crds_get_crd(crds, 0, FX1);
  float *fyy1 = ggcm_mhd_crds_get_crd(crds, 1, FX1);
  float *fzz1 = ggcm_mhd_crds_get_crd(crds, 2, FX1);
  float *fdx1 = ggcm_mhd_crds_get_crd(crds, 0, FD1);
  float *fdy1 = ggcm_mhd_crds_get_crd(crds, 1, FD1);
  float *fdz1 = ggcm_mhd_crds_get_crd(crds, 2, FD1);

  for (int i = -2; i < ldims[0] + 2; i++) {
    fxx1[i] = MRC_CRDX(mrc_crds, i);
  }
  for (int i = -2; i < ldims[1] + 2; i++) {
    fyy1[i] = MRC_CRDY(mrc_crds, i);
  }
  for (int i = -2; i < ldims[2] + 2; i++) {
    fzz1[i] = MRC_CRDZ(mrc_crds, i);
  }

  for (int i = -1; i < ldims[0] + 1; i++) {
    fdx1[i] = 1. / (.5 * (MRC_CRDX(mrc_crds, i+1) - MRC_CRDX(mrc_crds, i-1)));
  }
  for (int i = -1; i < ldims[1] + 1; i++) {
    fdy1[i] = 1. / (.5 * (MRC_CRDY(mrc_crds, i+1) - MRC_CRDY(mrc_crds, i-1)));
  }
  for (int i = -1; i < ldims[2] + 1; i++) {
    fdz1[i] = 1. / (.5 * (MRC_CRDZ(mrc_crds, i+1) - MRC_CRDZ(mrc_crds, i-1)));
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen subclass "mrc"

struct ggcm_mhd_crds_gen_ops ggcm_mhd_crds_gen_mrc_ops = {
  .name       = "mrc",
  .run        = ggcm_mhd_crds_gen_mrc_run,
};
