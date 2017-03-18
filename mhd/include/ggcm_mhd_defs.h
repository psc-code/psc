
#ifndef GGCM_MHD_DEFS
#define GGCM_MHD_DEFS

// ----------------------------------------------------------------------
// number of ghost points

#define BND (2)

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
  PSI,
  EE = UU,
};

// for primitive fields

enum {
  VX = 1,
  VY = 2,
  VZ = 3,
  PP = 4,
  CMSV = 5,
  N_PRIMITIVE = 8,
};

// ----------------------------------------------------------------------
// macros to ease field access

#define RR_(U, i,j,k, p)   M3(U, RR , i,j,k, p)
#define RVX_(U, i,j,k, p)  M3(U, RVX, i,j,k, p)
#define RVY_(U, i,j,k, p)  M3(U, RVY, i,j,k, p)
#define RVZ_(U, i,j,k, p)  M3(U, RVZ, i,j,k, p)
#define EE_(U, i,j,k, p)   M3(U, EE , i,j,k, p)
#define BX_(U, i,j,k, p)   M3(U, BX , i,j,k, p)
#define BY_(U, i,j,k, p)   M3(U, BY , i,j,k, p)
#define BZ_(U, i,j,k, p)   M3(U, BZ , i,j,k, p)
#define UU_(U, i,j,k, p)   M3(U, UU , i,j,k, p)

#define VX_(f, i,j,k, p)   M3(f, VX , i,j,k, p)
#define VY_(f, i,j,k, p)   M3(f, VY , i,j,k, p)
#define VZ_(f, i,j,k, p)   M3(f, VZ , i,j,k, p)
#define PP_(f, i,j,k, p)   M3(f, PP , i,j,k, p)

#define B0(b, d, i,j,k, p)  M3(b0, d, i,j,k, p)
#define B1(U, d, i,j,k, p)  M3(U, BX+d, i,j,k, p)
#define BT(U, d, i,j,k, p)  (b0 ? (B1(U, d, i,j,k, p) + B0(b0, d, i,j,k, p)) : B1(U, d, i,j,k, p))

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
  NR_CRDS,
};

// ----------------------------------------------------------------------
// assorted macros

#define sqr(x) ((x)*(x))

#endif
