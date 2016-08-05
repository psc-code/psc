
#include "ggcm_mhd_crds_gen_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_crds_private.h"

#include <mrc_domain.h>
#include <assert.h>
#include <string.h>

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen_mrc_run
//
// initializes the Fortran common block coord arrays FX1, FD1 from
// using C grid generation

static void
ggcm_mhd_crds_gen_mrc_run(struct ggcm_mhd_crds_gen *gen, struct ggcm_mhd_crds *crds)
{
}

// ----------------------------------------------------------------------
// ggcm_mhd_crds_gen subclass "mrc"

struct ggcm_mhd_crds_gen_ops ggcm_mhd_crds_gen_mrc_ops = {
  .name       = "mrc",
  .run        = ggcm_mhd_crds_gen_mrc_run,
};
