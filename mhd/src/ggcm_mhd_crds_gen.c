
#include "ggcm_mhd_crds_gen_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds_private.h"

#include <assert.h>

// ======================================================================
// ggcm_mhd_crds_gen class

#define ggcm_mhd_crds_gen_ops(gen) ((struct ggcm_mhd_crds_gen_ops *)((gen)->obj.ops))

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_run_aux_default
//
// default routine for calculating the remaining coordinate arrays

static void
ggcm_mhd_crds_gen_run_aux_default(struct ggcm_mhd_crds_gen *gen,
				  struct ggcm_mhd_crds *crds)
{
  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(crds->domain, 0, &info);
  int *im = info.ldims;

  float *fxx1 = ggcm_mhd_crds_get_crd(crds, 0, FX1);
  float *fyy1 = ggcm_mhd_crds_get_crd(crds, 1, FX1);
  float *fzz1 = ggcm_mhd_crds_get_crd(crds, 2, FX1);
  float *fxx2 = ggcm_mhd_crds_get_crd(crds, 0, FX2);
  float *fyy2 = ggcm_mhd_crds_get_crd(crds, 1, FX2);
  float *fzz2 = ggcm_mhd_crds_get_crd(crds, 2, FX2);

  for (int i = -2; i < im[0] + 2; i++) {
    fxx2[i] = sqr(fxx1[i]);
  }
  for (int i = -2; i < im[1] + 2; i++) {
    fyy2[i] = sqr(fyy1[i]);
  }
  for (int i = -2; i < im[2] + 2; i++) {
    fzz2[i] = sqr(fzz1[i]);
  }

  float *bdx1 = ggcm_mhd_crds_get_crd(crds, 0, BD1);
  float *bdy1 = ggcm_mhd_crds_get_crd(crds, 1, BD1);
  float *bdz1 = ggcm_mhd_crds_get_crd(crds, 2, BD1);
  float *bdx2 = ggcm_mhd_crds_get_crd(crds, 0, BD2);
  float *bdy2 = ggcm_mhd_crds_get_crd(crds, 1, BD2);
  float *bdz2 = ggcm_mhd_crds_get_crd(crds, 2, BD2);
  float *bdx3 = ggcm_mhd_crds_get_crd(crds, 0, BD3);
  float *bdy3 = ggcm_mhd_crds_get_crd(crds, 1, BD3);
  float *bdz3 = ggcm_mhd_crds_get_crd(crds, 2, BD3);
  float *bdx4 = ggcm_mhd_crds_get_crd(crds, 0, BD4);
  float *bdy4 = ggcm_mhd_crds_get_crd(crds, 1, BD4);
  float *bdz4 = ggcm_mhd_crds_get_crd(crds, 2, BD4);

  for (int i = -2; i < im[0] + 1; i++) {
    bdx1[i] = 1.f / (fxx1[i+1] - fxx1[i]);
    bdx4[i] = 1.f / (fxx1[i+1] - fxx1[i]);
  }
  for (int i = -1; i < im[0] + 1; i++) {
    bdx2[i] = .5f * (fxx1[i+1] - fxx1[i-1]);
    bdx3[i] = 1.f / bdx2[i];
  }
      
  for (int i = -2; i < im[1] + 1; i++) {
    bdy1[i] = 1.f / (fyy1[i+1] - fyy1[i]);
    bdy4[i] = 1.f / (fyy1[i+1] - fyy1[i]);
  }
  for (int i = -1; i < im[1] + 1; i++) {
    bdy2[i] = .5f * (fyy1[i+1] - fyy1[i-1]);
    bdy3[i] = 1.f / bdy2[i];
  }
#if 0
  { int i = -2;
    bdy2[i] = fyy1[i+1] - fyy1[i]; }
  { int i = im[1] + 1;
    bdy2[i] = fyy1[i] - fyy1[i-1]; }
#endif

  for (int i = -2; i < im[2] + 1; i++) {
    bdz1[i] = 1.f / (fzz1[i+1] - fzz1[i]);
    bdz4[i] = 1.f / (fzz1[i+1] - fzz1[i]);
  }
  for (int i = -1; i < im[2] + 1; i++) {
    bdz2[i] = .5f * (fzz1[i+1] - fzz1[i-1]);
    bdz3[i] = 1.f / bdz2[i];
  }
#if 0
  { int i = -2;
    bdz2[i] = fzz1[i+1] - fzz1[i]; }
  { int i = im[2] + 1;
    bdz2[i] = fzz1[i] - fzz1[i-1]; }
#endif
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_run

void
ggcm_mhd_crds_gen_run(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds)
{
  struct ggcm_mhd_crds_gen_ops *ops = ggcm_mhd_crds_gen_ops(gen);
  assert(ops && ops->run);

  ops->run(gen, crds);
  if (ops && ops->run_aux) {
    ops->run_aux(gen, crds);
  } else {
    ggcm_mhd_crds_gen_run_aux_default(gen, crds);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_init

static void
ggcm_mhd_crds_gen_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_crds_gen, &ggcm_mhd_crds_gen_c_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_crds_gen, &ggcm_mhd_crds_gen_mrc_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen class

struct mrc_class_ggcm_mhd_crds_gen mrc_class_ggcm_mhd_crds_gen = {
  .name             = "ggcm_mhd_crds_gen",
  .size             = sizeof(struct ggcm_mhd_crds_gen),
  .init             = ggcm_mhd_crds_gen_init,
};

