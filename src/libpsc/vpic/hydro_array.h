
#ifndef HYDRO_ARRAY_H
#define HYDRO_ARRAY_H

#include "grid.h"

// ======================================================================
// HydroArray

struct HydroArray : hydro_array_t {
  HydroArray(Grid g);
  ~HydroArray();
};

// ----------------------------------------------------------------------
// copied from hydro_array.c, converted from new/delete -> ctor

inline void
hydro_array_ctor(hydro_array_t * ha, grid_t * g ) {
  if( !g ) ERROR(( "NULL grid" ));
  MALLOC_ALIGNED( ha->h, g->nv, 128 );
  ha->g = g;
  clear_hydro_array( ha );
}

inline void
hydro_array_dtor( hydro_array_t * ha ) {
  if( !ha ) return;
  FREE_ALIGNED( ha->h );
}

// ----------------------------------------------------------------------
// HydroArray implementation

inline HydroArray::HydroArray(Grid grid)
{
  hydro_array_ctor(this, grid.g_);
}

inline HydroArray::~HydroArray()
{
  hydro_array_dtor(this);
}


#endif

