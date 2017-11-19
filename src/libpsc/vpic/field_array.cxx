
#include "field_array.h"

#include <mrc_common.h>

inline float FieldArray::operator()(int m, int i, int j, int k) const
{
  float *f_ = &f[0].ex;
  return f_[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
}

inline float& FieldArray::operator()(int m, int i, int j, int k)
{
  float *f_ = &f[0].ex;
  return f_[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
}

struct Field3D {
  Field3D(FieldArray& fa)
    : sx_(fa.g->nx + 2), sy_(fa.g->ny + 2),
      f_(fa.data())
  {
  }

  int voxel(int i, int j, int k) const
  {
    return i + sx_ * (j + sy_ * (k));
  }

  field_t& operator()(FieldArray &fa, int i, int j, int k)
  {
    return fa.f[voxel(i,j,k)];
  }
  
  field_t operator()(FieldArray &fa, int i, int j, int k) const
  {
    return fa.f[voxel(i,j,k)];
  }
  
  float& operator()(int m, int i, int j, int k)
  {
    return f_[m + FieldArray::N_COMP * voxel(i,j,k)];
  }
  
  float operator()(int m, int i, int j, int k) const
  {
    return f_[m + FieldArray::N_COMP * voxel(i,j,k)];
  }
  
  int sx_, sy_;
  float * RESTRICT f_;
};

#define EX(i,j,k)  F(EX , i,j,k)
#define EY(i,j,k)  F(EY , i,j,k)
#define EZ(i,j,k)  F(EZ , i,j,k)
#define CBX(i,j,k) F(CBX, i,j,k)
#define CBY(i,j,k) F(CBY, i,j,k)
#define CBZ(i,j,k) F(CBZ, i,j,k)

#define DECLARE_STENCIL()						\
  Field3D F(*this);							\
  const int   nx   = g->nx;						\
  const int   ny   = g->ny;						\
  const int   nz   = g->nz;						\
  									\
  const float px   = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;		\
  const float py   = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;		\
  const float pz   = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0

#define f(i,j,k) f[ VOXEL(i,j,k, nx,ny,nz) ]

#define UPDATE_CBX() CBX(i,j,k) -= (py*(EZ(i,j+1,k) - EZ(i,j,k)) - pz*(EY(i,j,k+1) - EY(i,j,k)))
#define UPDATE_CBY() CBY(i,j,k) -= (pz*(EX(i,j,k+1) - EX(i,j,k)) - px*(EZ(i+1,j,k) - EZ(i,j,k)))
#define UPDATE_CBZ() CBZ(i,j,k) -= (px*(EY(i+1,j,k) - EY(i,j,k)) - py*(EX(i,j+1,k) - EX(i,j,k)))

void FieldArray::advanceB_interior(float frac)
{
  DECLARE_STENCIL();

  for (int k = 1; k <= nz; k++) {
    for (int j = 1; j <= ny; j++) {
      for (int i = 1; i <= nx; i++) {
	UPDATE_CBX(); UPDATE_CBY(); UPDATE_CBZ();
      }
    }
  }
}

void FieldArray::advanceB(float frac)
{
  advanceB_interior(frac);

  DECLARE_STENCIL();
  
  // leftover bx
  { int i = nx + 1;
    for (int k = 1; k <= nz; k++) {
      for (int j = 1; j <= ny; j++) {
	UPDATE_CBX();
      }
    }
  }

  // leftover by
  { int j = ny + 1;
    for (int k = 1; k <= nz; k++) {
      for (int i = 1; i <= nx; i++) {
	UPDATE_CBY();
      }
    }
  }

  // leftover bz
  { int k = nz + 1;
    for (int j = 1; j <= ny; j++) {
      for (int i = 1; i <= nx; i++) {
	UPDATE_CBZ();
      }
    }
  }
  
  local_adjust_norm_b( f, g );
}

