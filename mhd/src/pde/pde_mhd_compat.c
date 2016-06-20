
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

#define FX1X(ix) PDE_CRDX_CC(ix)
#define FX1Y(iy) PDE_CRDY_CC(iy)
#define FX1Z(iz) PDE_CRDZ_CC(iz)

// FD1 is really different if 'legacy_fd1' is used
#if 1
#define FD1X(ix) (s_fd1x[ix])
#define FD1Y(iy) (s_fd1y[iy])
#define FD1Z(iz) (s_fd1z[iz])
#else
#define FD1X(ix) PDE_INV_DX(ix)
#define FD1Y(iy) PDE_INV_DY(iy)
#define FD1Z(iz) PDE_INV_DZ(iz)
#endif

#define BD1X(ix) PDE_INV_DXF(ix+1)
#define BD1Y(iy) PDE_INV_DYF(iy+1)
#define BD1Z(iz) PDE_INV_DZF(iz+1)

#define BD2X(ix) PDE_DX(ix)
#define BD2Y(iy) PDE_DY(iy)
#define BD2Z(iz) PDE_DZ(iz)

#define BD3X(ix) PDE_INV_DX(ix)
#define BD3Y(iy) PDE_INV_DY(iy)
#define BD3Z(iz) PDE_INV_DZ(iz)

#define BD4X(ix) BD1X(ix)
#define BD4Y(iy) BD1Y(iy)
#define BD4Z(iz) BD1Z(iz)

#endif
