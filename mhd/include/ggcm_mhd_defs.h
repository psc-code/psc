
#ifndef GGCM_MHD_DEFS
#define GGCM_MHD_DEFS

// ----------------------------------------------------------------------
// number of ghost points

#define BND (2)

// ----------------------------------------------------------------------
// fields
// FIXME, these are dependent on what openggcm does
// indices based on mhd-corea.for

enum {
  _RR1,
  _RV1X,
  _RV1Y,
  _RV1Z,
  _UU1,
  __B1X,
  __B1Y,
  __B1Z,

  __RR2,
  __RV2X,
  __RV2Y,
  __RV2Z,
  __UU2,
  __B2X,
  __B2Y,
  __B2Z,

  __YMASK, // 16
  __ZMASK,
  __CMSV,
  
  __RR, // 19
  __PP,
  __VX,
  __VY,
  __VZ,
  __BX,
  __BY,
  __BZ,

  __TMP1, // 27
  __TMP2,
  __TMP3,
  __TMP4,

  __FLX, // 31
  __FLY,
  __FLZ,

  __CX, // 34
  __CY,
  __CZ,
  __XTRA1, // 37
  __XTRA2,
  __RESIS,
  __CURRX, // 40
  __CURRY,
  __CURRZ,
  __RMASK, // 43
  _NR_FLDS,
};

// for the state vector

enum {
  RR,
  RVX,
  RVY,
  RVZ,
  UU,
  BX,
  BY,
  BZ,
  EE = UU,
};

// for primitive fields

enum {
  VX = 1,
  VY = 2,
  VZ = 3,
  PP = 4,
  CMSV = 5,
};

// ----------------------------------------------------------------------
// macros to ease field access

#define BXYZ(f,m, ix,iy,iz) F3(f, BX+(m), ix,iy,iz)

#define RR(U, i,j,k)   F3(U, RR , i,j,k)
#define RVX(U, i,j,k)  F3(U, RVX, i,j,k)
#define RVY(U, i,j,k)  F3(U, RVY, i,j,k)
#define RVZ(U, i,j,k)  F3(U, RVZ, i,j,k)
#define EE(U, i,j,k)   F3(U, EE , i,j,k)
#define BX(U, i,j,k)   F3(U, BX , i,j,k)
#define BY(U, i,j,k)   F3(U, BY , i,j,k)
#define BZ(U, i,j,k)   F3(U, BZ , i,j,k)
#define UU(U, i,j,k)   F3(U, UU , i,j,k)

#define VX(f, ix,iy,iz) F3(f, VX, ix,iy,iz)
#define VY(f, ix,iy,iz) F3(f, VY, ix,iy,iz)
#define VZ(f, ix,iy,iz) F3(f, VZ, ix,iy,iz)
#define PP(f, ix,iy,iz) F3(f, PP, ix,iy,iz)

// ----------------------------------------------------------------------
// coordinates

enum {
  FX1, // x_{i+1/2} cell-centered, [-1:mx-1]
  FX2, // same as FX1, squared
  FD1, // 1 / (x_{i+1} - x_{i}) = 1 / (.5*(FX1[i+1] - FX1[i-1]))
  BD1, // 1 / (x_{i+3/2} - x_{i+1/2}) = 1 / (FX1[i+1] - FX1[i])
  BD2, // x_{i+1} - x_{i} = .5*(FX1[i+1] - FX1[i-1])
  BD3, // 1 / BD2 == FD1
  BD4, // == BD1
  NR_CRDS, // FIXME, update from Fortran
};

// ----------------------------------------------------------------------
// assorted macros

#define sqr(x) ((x)*(x))

#endif
