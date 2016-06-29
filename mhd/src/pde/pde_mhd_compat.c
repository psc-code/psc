
#ifndef PDE_MHD_COMPAT_C
#define PDE_MHD_COMPAT_C

#include "ggcm_mhd_crds.h"

#include "pde/pde_mhd_setup.c"

// ======================================================================
// compatibility functionality for original openggcm
//
// - support FD1 coordinate, which should just be the (inverse) cell size,
//   but is calculated from the derivative of the coordinate mapping function,
//   and hence not exactly the same

static mrc_fld_data_t s_mhd_time;

static float *s_fd1x, *s_fd1y, *s_fd1z;

// we keep a static copy of the mhd->mrc_fld (as fld3d_t) around,
// in order to have fld3d_make_tmp() actually provide the temporaries in their
// original location. also used for fortran interfacing
static fld3d_t s_p_f;

// FD1 is really different if 'legacy_fd1' is used
#if 1
#define FD1X(i) (s_fd1x[i])
#define FD1Y(j) (s_fd1y[j])
#define FD1Z(k) (s_fd1z[k])
#else
#define FD1X(i) PDE_INV_DX(i)
#define FD1Y(j) PDE_INV_DY(j)
#define FD1Z(k) PDE_INV_DZ(k)
#endif

static void
pde_mhd_compat_setup(struct ggcm_mhd *mhd)
{
  s_fd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, FD1);
  s_fd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, FD1);
  s_fd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, FD1);

  // for temp / Fortran interface, in which case we only have one patch
  fld3d_setup(&s_p_f, mhd->fld);
  fld3d_get(&s_p_f, 0);
}

#if OPT_TMP == OPT_TMP_COMPAT

#define fld3d_setup_tmp_compat(p_tmp, n_comps, m) do {	\
    if (!(p_tmp)->arr_off) {				\
      fld3d_setup_view(p_tmp, s_p_f, m);		\
    }							\
  } while (0)

#else

#define fld3d_setup_tmp_compat(p_tmp, n_comps, m) do {	\
    if (!(p_tmp)->arr_off) {				\
      fld3d_setup_tmp(p_tmp, n_comps);			\
    }							\
  } while (0)

#endif

#endif
