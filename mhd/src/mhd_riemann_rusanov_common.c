
#include "mhd_riemann_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"

#include <math.h>

#include "pde/pde_setup.c"
#include "pde/pde_mhd_riemann.c"

// ----------------------------------------------------------------------
// mhd_riemann_rusanov_run

static void
mhd_riemann_rusanov_run(struct mhd_riemann *riemann, struct mrc_fld *F,
			struct mrc_fld *U_l, struct mrc_fld *U_r,
			struct mrc_fld *W_l, struct mrc_fld *W_r,
			int ldim, int l, int r, int dim)
{
  Gamma = riemann->mhd->par.gamm;

  for (int i = -l; i < ldim + r; i++) {
    fluxes_rusanov_sc(&F1(F, 0, i), &F1(U_l, 0, i), &F1(U_r, 0, i),
		      &F1(W_l, 0, i), &F1(W_r, 0, i));
  }
}

