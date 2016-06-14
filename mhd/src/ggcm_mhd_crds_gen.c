
#include "ggcm_mhd_crds_gen_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds_private.h"

#include <assert.h>

// ======================================================================
// ggcm_mhd_crds_gen class

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_run
//
// initializes the ggcm_mhd_crds arrays form mrc_crds

void
ggcm_mhd_crds_gen_run(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds)
{
  struct mrc_crds *mrc_crds = mrc_domain_get_crds(crds->domain);

  for (int p = 0; p < mrc_domain_nr_patches(crds->domain); p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(crds->domain, p, &info);
    int *ldims = info.ldims;
    for (int d = 0; d < 3; d++) {
      // yes, these are labeled with 'x', but we loop over all dimensions here
      float *fxx1 = ggcm_mhd_crds_get_crd_p(crds, d, FX1, p);
      float *fdx1 = ggcm_mhd_crds_get_crd_p(crds, d, FD1, p);
      float *fxx2 = ggcm_mhd_crds_get_crd_p(crds, d, FX2, p);
      float *bdx1 = ggcm_mhd_crds_get_crd_p(crds, d, BD1, p);
      float *bdx2 = ggcm_mhd_crds_get_crd_p(crds, d, BD2, p);
      float *bdx3 = ggcm_mhd_crds_get_crd_p(crds, d, BD3, p);
      float *bdx4 = ggcm_mhd_crds_get_crd_p(crds, d, BD4, p);
      int sw = mrc_crds->sw;
      
      for (int i = -sw; i < ldims[d] + sw; i++) {
	fxx1[i] = MRC_MCRD(mrc_crds, d, i, p);
      }
      
      // have to move one in on both sides
      for (int i = -sw + 1; i < ldims[d] + sw - 1; i++) {
	if (gen->legacy_fd1) {
	  int off = info.off[d];
	  fdx1[i] = 1.f / MRC_D2(mrc_crds->global_crd[d], i + off, 1);
	} else {
	  fdx1[i] = 1.f / (.5f * (MRC_MCRD(mrc_crds, d, i+1, p) - MRC_MCRD(mrc_crds, d, i-1, p)));
	}
      }

      for (int i = -sw; i < ldims[d] + sw; i++) {
	fxx2[i] = sqr(MRC_CRD(mrc_crds, d, i));
      }

      for (int i = -sw; i < ldims[d] + sw - 1; i++) {
	bdx1[i] = 1.f / (MRC_CRD(mrc_crds, d, i+1) - MRC_CRD(mrc_crds, d, i));
	bdx4[i] = 1.f / (MRC_CRD(mrc_crds, d, i+1) - MRC_CRD(mrc_crds, d, i));
      }

      for (int i = -sw; i < ldims[d] + sw; i++) {
	float old = bdx2[i];
	if (i == -sw || i == ldims[d] + sw - 1) {
	  // we could just use this always, but it gives finite precision
	  // errors compared to the version below, which is how fortran did it
	  bdx2[i] = MRC_MCRD_NC(mrc_crds, d, i+1, 0) - MRC_MCRD_NC(mrc_crds, d, i, 0);
	} else {
	  bdx2[i] = .5f * (MRC_MCRD(mrc_crds, d, i+1, 0) - MRC_MCRD(mrc_crds, d, i-1, 0));
	}
	bdx3[i] = 1.f / bdx2[i];
      }
    }
  }

  if (strcmp(mrc_domain_type(crds->domain), "simple") == 0) {
    for (int d = 0; d < 3; d++) {
      struct mrc_fld *global_x = crds->global_f1[d];
      mrc_f1_foreach(global_x, i, 0, 0) {
	MRC_F1(global_x, 0, i) = MRC_D2(mrc_crds->global_crd[d], i, 0);
      } mrc_f1_foreach_end;
    }
  }
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
};

