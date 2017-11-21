
#ifndef VPIC_INTERPOLATOR_H
#define VPIC_INTERPOLATOR_H

#include "grid.h"

// ======================================================================
// VpicInterpolator

struct VpicInterpolator : interpolator_array_t {
  VpicInterpolator(Grid g);
  ~VpicInterpolator();
};

// ----------------------------------------------------------------------
// copied from interpolator_array.c, converted from new/delete -> ctor

inline void
interpolator_array_ctor(interpolator_array_t * ia, grid_t * g ) {
  if( !g ) ERROR(( "NULL grid" ));
  MALLOC_ALIGNED( ia->i, g->nv, 128 );
  CLEAR( ia->i, g->nv );
  ia->g = g;
}

inline void
interpolator_array_dtor( interpolator_array_t * ia ) {
  if( !ia ) return;
  FREE_ALIGNED( ia->i );
}

// ----------------------------------------------------------------------
// VpicInterpolator implementation

inline VpicInterpolator::VpicInterpolator(Grid grid)
{
  interpolator_array_ctor(this, grid.getGrid_t());
}

inline VpicInterpolator::~VpicInterpolator()
{
  interpolator_array_dtor(this);
}



#endif

