
#ifndef INTERPOLATOR_ARRAY_H
#define INTERPOLATOR_ARRAY_H

#include "grid.h"

// ======================================================================
// InterpolatorArray

struct InterpolatorArray : interpolator_array_t {
  InterpolatorArray(Grid g);
  ~InterpolatorArray();
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
// InterpolatorArray implementation

inline InterpolatorArray::InterpolatorArray(Grid grid)
{
  interpolator_array_ctor(this, grid.getGrid_t());
}

inline InterpolatorArray::~InterpolatorArray()
{
  interpolator_array_dtor(this);
}



#endif

