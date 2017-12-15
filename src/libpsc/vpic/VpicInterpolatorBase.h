
#ifndef VPIC_INTERPOLATOR_BASE_H
#define VPIC_INTERPOLATOR_BASE_H

// ======================================================================
// VpicInterpolatorBase

template<class G>
struct VpicInterpolatorBase : interpolator_array_t {
  typedef G Grid;
  typedef interpolator_t Element;

  enum {
    EX        = 0,
    DEXDY     = 1,
    DEXDZ     = 2,
    D2DEXDYDZ = 3,
    EY        = 4,
    DEYDZ     = 5,
    DEYDX     = 6,
    D2DEYDZDX = 7,
    EZ        = 8,
    DEZDX     = 9,
    DEZDY     = 10,
    D2DEZDXDY = 11,
    CBX       = 12,
    DCBXDX    = 13,
    CBY       = 14,
    DCBYDY    = 15,
    CBZ       = 16,
    DCBZDZ    = 17,
    
    N_COMP = sizeof(interpolator_t) / sizeof(float),
  };

  static VpicInterpolatorBase* create(Grid *grid)
  {
    return static_cast<VpicInterpolatorBase*>(new_interpolator_array(grid));
  }
  
  static void destroy(VpicInterpolatorBase* interpolator)
  {
    delete_interpolator_array(interpolator);
  }
  
  Element operator[](int idx) const
  {
    return i[idx];
  }

  Element& operator[](int idx)
  {
    return i[idx];
  }

  Element* data()
  {
    return i;
  }
};

#if 0
// ----------------------------------------------------------------------
// copied from interpolator_array.c, converted from new/delete -> ctor

inline void
interpolator_array_ctor(interpolator_array_t * ia, grid_t * g) {
  if( !g ) ERROR(( "NULL grid" ));
  MALLOC_ALIGNED( ia->i, g->nv, 128 );
  CLEAR( ia->i, g->nv );
  ia->g = g;
}

inline void
interpolator_array_dtor(interpolator_array_t * ia) {
  if( !ia ) return;
  FREE_ALIGNED( ia->i );
}
#endif

#endif

