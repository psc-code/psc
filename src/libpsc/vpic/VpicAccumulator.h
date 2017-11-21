
#ifndef VPIC_ACCUMULATOR_H
#define VPIC_ACCUMULATOR_H

#include "grid.h"

// ======================================================================
// VpicAccumulator

struct VpicAccumulator : accumulator_array_t {
  VpicAccumulator(Grid g);
  ~VpicAccumulator();
};

// ----------------------------------------------------------------------
// copied from accumulator_array.c, converted from new/delete -> ctor

static int
aa_n_pipeline(void) {
  int                       n = serial.n_pipeline;
  if( n<thread.n_pipeline ) n = thread.n_pipeline;
  return n; /* max( {serial,thread,spu}.n_pipeline ) */
}

inline void
accumulator_array_ctor(accumulator_array_t * aa, grid_t * g ) {
  if( !g ) ERROR(( "Bad grid."));
  aa->n_pipeline = aa_n_pipeline();
  aa->stride     = POW2_CEIL(g->nv,2);
  aa->g          = g;
  MALLOC_ALIGNED( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride, 128 );
  CLEAR( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride );
}

inline void
accumulator_array_dtor( accumulator_array_t * aa ) {
  if( !aa ) return;
  FREE_ALIGNED( aa->a );
}

// ----------------------------------------------------------------------
// VpicAccumulator implementation

inline VpicAccumulator::VpicAccumulator(Grid grid)
{
  accumulator_array_ctor(this, grid.getGrid_t());
}

inline VpicAccumulator::~VpicAccumulator()
{
  accumulator_array_dtor(this);
}



#endif

