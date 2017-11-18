
#include "field_array.h"

#include <mrc_common.h>

#define DECLARE_STENCIL()                                       \
  const int   nx   = g->nx;                                     \
  const int   ny   = g->ny;                                     \
  const int   nz   = g->nz;                                     \
                                                                \
  const float px   = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;    \
  const float py   = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;    \
  const float pz   = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0

#define f(i,j,k) f[ VOXEL(i,j,k, nx,ny,nz) ]

#define UPDATE_CBX() f(i,j,k).cbx -= (py*(f(i,j+1,k).ez - f(i,j,k).ez) - pz*(f(i,j,k+1).ey - f(i,j,k).ey))
#define UPDATE_CBY() f(i,j,k).cby -= (pz*(f(i,j,k+1).ex - f(i,j,k).ex) - px*(f(i+1,j,k).ez - f(i,j,k).ez))
#define UPDATE_CBZ() f(i,j,k).cbz -= (px*(f(i+1,j,k).ey - f(i,j,k).ey) - py*(f(i,j+1,k).ex - f(i,j,k).ex))

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

