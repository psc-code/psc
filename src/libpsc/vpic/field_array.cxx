
#include "field_array.h"

#include <mrc_common.h>

#define DECLARE_STENCIL()                                       \
  /**/  field_t * ALIGNED(128) f = fa->f;			\
  const grid_t  *              g = fa->g;			\
                                                                \
  const int   nx   = g->nx;                                     \
  const int   ny   = g->ny;                                     \
  const int   nz   = g->nz;                                     \
                                                                \
  const float px   = (nx>1) ? frac*g->cvac*g->dt*g->rdx : 0;    \
  const float py   = (ny>1) ? frac*g->cvac*g->dt*g->rdy : 0;    \
  const float pz   = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

#define UPDATE_CBX() f(x,y,z).cbx -= (py*(f(x,y+1,z).ez - f(x,y,z).ez) - pz*(f(x,y,z+1).ey - f(x,y,z).ey))
#define UPDATE_CBY() f(x,y,z).cby -= (pz*(f(x,y,z+1).ex - f(x,y,z).ex) - px*(f(x+1,y,z).ez - f(x,y,z).ez))
#define UPDATE_CBZ() f(x,y,z).cbz -= (px*(f(x+1,y,z).ey - f(x,y,z).ey) - py*(f(x,y+1,z).ex - f(x,y,z).ex))

static void
advance_b_interior(field_array_t *fa, float frac)
{
  DECLARE_STENCIL();

  for (int z = 1; z <= nz; z++) {
    for (int y = 1; y <= ny; y++) {
      for (int x = 1; x <= nx; x++) {
	UPDATE_CBX(); UPDATE_CBY(); UPDATE_CBZ();
      }
    }
  }
}


static void
advance_b_(field_array_t* fa, float frac)
{
  advance_b_interior(fa, frac);

  DECLARE_STENCIL();
  int x, y, z;
  
  // Do left over bx
  x = nx+1;
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      UPDATE_CBX();
    }
  }

  // Do left over by
  y = ny+1;
  for( z=1; z<=nz; z++ ) {
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBY();
    }
  }

  // Do left over bz
  z = nz+1;
  for( y=1; y<=ny; y++ ) {
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBZ();
    }
  }

  local_adjust_norm_b( f, g );
}

void FieldArray::advanceB(double frac)
{
  advance_b_(this, frac);
}

