
#ifndef GGCM_MHD_DEFS
#define GGCM_MHD_DEFS

// ----------------------------------------------------------------------
// number of ghost points

#define BND (2)

// ----------------------------------------------------------------------
// fields
// FIXME, these are dependent on what openggcm does

enum {
  _RR1,
  _RV1X,
  _RV1Y,
  _RV1Z,
  _UU1,
  _B1X,
  _B1Y,
  _B1Z,
};

// ----------------------------------------------------------------------
// macros to ease field access

#if 0
#define B1XYZ(f,m, ix,iy,iz) MRC_F3(f, _B1X+(m), ix,iy,iz)
#else // GGCM staggering
#define B1XYZ(f,m, ix,iy,iz) MRC_F3(f, _B1X+(m),		\
				    (ix) - ((m) == 0),		\
				    (iy) - ((m) == 1),		\
				    (iz) - ((m) == 2))
#endif

#define B1X(f, ix,iy,iz) B1XYZ(f, 0, ix,iy,iz)
#define B1Y(f, ix,iy,iz) B1XYZ(f, 1, ix,iy,iz)
#define B1Z(f, ix,iy,iz) B1XYZ(f, 2, ix,iy,iz)

// ----------------------------------------------------------------------
// coordinates

enum {
  FX1, // node-centered, [-1:mx-1]
  FX2, // same as FX1, squared
  FD1,
  BD1,
  BD2,
  BD3,
  BD4,
  NR_CRDS, // FIXME, update from Fortran
};

// ----------------------------------------------------------------------
// assorted macros

#define sqr(x) ((x)*(x))

#endif
