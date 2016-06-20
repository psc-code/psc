
#ifndef PDE_MHD_COMPAT_C
#define PDE_MHD_COMPAT_C

// ======================================================================
// compatibility functionality for original openggcm
//
// - support FD1 coordinate, which should just be the (inverse) cell size,
//   but is calculated from the derivative of the coordinate mapping function,
//   and hence not exactly the same

static mrc_fld_data_t s_mhd_time;

static float *s_fd1x, *s_fd1y, *s_fd1z;

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

#endif
