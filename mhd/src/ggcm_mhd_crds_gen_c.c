
#include "ggcm_mhd_crds_gen_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds_private.h"
#include "mrc_crds_gen.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_c_run
//
// initializes the Fortran common block coord arrays FX1, FD1 from
// using C grid generation

static void
ggcm_mhd_crds_gen_c_run(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds)
{
  MPI_Comm comm = ggcm_mhd_crds_gen_comm(gen);

  int gdims[3];
  mrc_domain_get_global_dims(crds->domain, gdims);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(crds->domain, 0, &info);
  int *ldims = info.ldims;

  float xl[3], xh[3];
  struct mrc_crds *mrc_crds = mrc_domain_get_crds(crds->domain);
  mrc_crds_get_param_float3(mrc_crds, "l", xl);
  mrc_crds_get_param_float3(mrc_crds, "h", xh);

  for (int d = 0; d < 3; d++) {
    double *xx = malloc((gdims[d] + 2*BND + 1) * sizeof(double));
    double *dx = malloc((gdims[d] + 2*BND + 1) * sizeof(double));

    struct mrc_crds_gen *gen = mrc_crds_gen_create(comm);
    if (d == 0) {
      mrc_crds_gen_set_type(gen, "ggcm_x");
    } else {
      mrc_crds_gen_set_type(gen, "ggcm_yz");
    }
    char name[10];
    sprintf(name, "grid%c", 'x' + d);
    mrc_crds_gen_set_name(gen, name);
    mrc_crds_gen_set_param_obj(gen, "crds", mrc_crds);
    mrc_crds_gen_set_param_int(gen, "d", d);
    mrc_crds_gen_set_from_options(gen);
    mrc_crds_gen_setup(gen);
    mrc_crds_gen_view(gen);
    mrc_crds_gen_run(gen, xx + BND, dx + BND);
    mrc_crds_gen_destroy(gen);

    // shift to beginning of local domain
    double *xxl = xx + info.off[d] + 2;

    // copy local part into fortran arrays
    float *fxx1 = ggcm_mhd_crds_get_crd(crds, d, FX1);
    float *fdx1 = ggcm_mhd_crds_get_crd(crds, d, FD1);

    for (int i = -2; i < ldims[d] + 2; i++) {
      fxx1[i] = xxl[i];
    }

    for (int i = -1; i < ldims[d] + 1; i++) {
      fdx1[i] = 1. / (.5 * (xxl[i + 1] - xxl[i - 1]));
    }

    free(xx);
    free(dx);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen subclass "c"

struct ggcm_mhd_crds_gen_ops ggcm_mhd_crds_gen_c_ops = {
  .name       = "c",
  .run        = ggcm_mhd_crds_gen_c_run,
};
