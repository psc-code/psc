
#include "field_array.h"

#include <mrc_common.h>

#if 0
#define EX(i,j,k) (*this)(EX, i,j,k)
#define EY(i,j,k) (*this)(EY, i,j,k)
#define EZ(i,j,k) (*this)(EZ, i,j,k)
#define CBX(i,j,k) (*this)(CBX, i,j,k)
#define CBY(i,j,k) (*this)(CBY, i,j,k)
#define CBZ(i,j,k) (*this)(CBZ, i,j,k)
#endif

#if 0
#define EX(i,j,k) ex(i,j,k)
#define EY(i,j,k) ey(i,j,k)
#define EZ(i,j,k) ez(i,j,k)
#define CBX(i,j,k) cbx(i,j,k)
#define CBY(i,j,k) cby(i,j,k)
#define CBZ(i,j,k) cbz(i,j,k)
#endif

#if 1
#define EX(i,j,k) f(i,j,k).ex
#define EY(i,j,k) f(i,j,k).ey
#define EZ(i,j,k) f(i,j,k).ez
#define CBX(i,j,k) f(i,j,k).cbx
#define CBY(i,j,k) f(i,j,k).cby
#define CBZ(i,j,k) f(i,j,k).cbz
#endif

#define VOXEL(x,y,z, nx,ny,nz) ((x) + ((nx)+2)*((y) + ((ny)+2)*(z)))
inline float FieldArray::operator()(int m, int i, int j, int k) const
{
  float *ff = &f[0].ex;
  return ff[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
}

inline float& FieldArray::operator()(int m, int i, int j, int k)
{
  float *ff = &f[0].ex;
  return ff[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
}

#define DECLARE_STENCIL()                                       \
  const int   nx   = g->nx;                                     \
  const int   ny   = g->ny;                                     \
  const int   nz   = g->nz;                                     \
                                                                \
  const float px   = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;    \
  const float py   = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;    \
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

