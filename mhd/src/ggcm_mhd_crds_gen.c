
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
  struct mrc_crds *mrc_crds = mrc_domain_get_crds(crds->domain);
  mrc_domain_get_local_patch_info(crds->domain, 0, &info);
  int *im = info.ldims;
  int sw = mrc_crds->sw;

  // FIXME: this is not bounds checked at all since sw is set from mrc_crds,
  //        NOT ggcm_mhd_crds, which is on the LHS OR RHS here
  for (int p = 0; p < mrc_domain_nr_patches(crds->domain); p++) {
    for (int d = 0; d < 3; d++) {
      // yes, these are labled with 'x', but we loops over all dimensions here
      float *fxx1 = ggcm_mhd_crds_get_crd_p(crds, d, FX1, p);
      float *fxx2 = ggcm_mhd_crds_get_crd_p(crds, d, FX2, p);
      float *bdx1 = ggcm_mhd_crds_get_crd_p(crds, d, BD1, p);
      float *bdx2 = ggcm_mhd_crds_get_crd_p(crds, d, BD2, p);
      float *bdx3 = ggcm_mhd_crds_get_crd_p(crds, d, BD3, p);
      float *bdx4 = ggcm_mhd_crds_get_crd_p(crds, d, BD4, p);
      
      for (int i = -sw; i < im[d] + sw; i++) {
	fxx2[i] = sqr(fxx1[i]);
      }
      for (int i = -sw; i < im[d] + sw - 1; i++) {
	bdx1[i] = 1.0 / (fxx1[i+1] - fxx1[i]);
	bdx4[i] = 1.0 / (fxx1[i+1] - fxx1[i]);
      }
      for (int i = -sw + 1; i < im[d] + sw - 1; i++) {
	bdx2[i] = 0.5 * (fxx1[i+1] - fxx1[i-1]);
	bdx3[i] = 1.0 / bdx2[i];
      }
    }
  }
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
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_crds_gen, &ggcm_mhd_crds_gen_mrc_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_crds_gen, x)
static struct param ggcm_mhd_crds_gen_descr[] = {
  { "legacy_fd1",      VAR(legacy_fd1),      PARAM_BOOL(0)     },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen class

struct mrc_class_ggcm_mhd_crds_gen mrc_class_ggcm_mhd_crds_gen = {
  .name             = "ggcm_mhd_crds_gen",
  .size             = sizeof(struct ggcm_mhd_crds_gen),
  .param_descr      = ggcm_mhd_crds_gen_descr,
  .init             = ggcm_mhd_crds_gen_init,
};

