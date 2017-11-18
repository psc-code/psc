
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
  const float pz   = (nz>1) ? frac*g->cvac*g->dt*g->rdz : 0;    \
                                                                \
  field_t * ALIGNED(16) f0;                                     \
  field_t * ALIGNED(16) fx, * ALIGNED(16) fy, * ALIGNED(16) fz; \
  int x, y, z

#define f(x,y,z) f[ VOXEL(x,y,z, nx,ny,nz) ]

#define INIT_STENCIL()  \
  f0 = &f(x,  y,  z  ); \
  fx = &f(x+1,y,  z  ); \
  fy = &f(x,  y+1,z  ); \
  fz = &f(x,  y,  z+1)

#define NEXT_STENCIL()                \
  f0++; fx++; fy++; fz++; x++;        \
  if( x>nx ) {                        \
    /**/       y++;            x = 1; \
    if( y>ny ) z++; if( y>ny ) y = 1; \
    INIT_STENCIL();                   \
  }
 
#define UPDATE_CBX() f0->cbx -= ( py*( fy->ez-f0->ez ) - pz*( fz->ey-f0->ey ) )
#define UPDATE_CBY() f0->cby -= ( pz*( fz->ex-f0->ex ) - px*( fx->ez-f0->ez ) )
#define UPDATE_CBZ() f0->cbz -= ( px*( fx->ey-f0->ey ) - py*( fy->ex-f0->ex ) )

static void
advance_b_pipeline(field_array_t *fa,
		   float frac,
                   int pipeline_rank,
                   int n_pipeline) {
  DECLARE_STENCIL();

  int n_voxel = nx * ny * nz;
  x = 1, y = 1, z = 1;

  INIT_STENCIL();
  for( ; n_voxel; n_voxel-- ) {
    UPDATE_CBX(); UPDATE_CBY(); UPDATE_CBZ();
    NEXT_STENCIL();
  }
}


static void
advance_b_(field_array_t* fa, float frac)
{
  advance_b_pipeline(fa, frac, 0, 1);

  // While the pipelines are busy, do surface fields

  DECLARE_STENCIL();
  
  // Do left over bx
  for( z=1; z<=nz; z++ ) {
    for( y=1; y<=ny; y++ ) {
      f0 = &f(nx+1,y,  z);
      fy = &f(nx+1,y+1,z);
      fz = &f(nx+1,y,  z+1);
      UPDATE_CBX();
    }
  }

  // Do left over by
  for( z=1; z<=nz; z++ ) {
    f0 = &f(1,ny+1,z);
    fx = &f(2,ny+1,z);
    fz = &f(1,ny+1,z+1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBY();
      f0++;
      fx++;
      fz++;
    }
  }

  // Do left over bz
  for( y=1; y<=ny; y++ ) {
    f0 = &f(1,y,  nz+1);
    fx = &f(2,y,  nz+1);
    fy = &f(1,y+1,nz+1);
    for( x=1; x<=nx; x++ ) {
      UPDATE_CBZ();
      f0++;
      fx++;
      fy++;
    }
  }

  local_adjust_norm_b( f, g );
}

void FieldArray::advanceB(double frac)
{
  advance_b_(this, frac);
}

