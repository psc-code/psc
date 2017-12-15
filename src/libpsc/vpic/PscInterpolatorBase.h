
#ifndef PSC_INTERPOLATOR_BASE_H
#define PSC_INTERPOLATOR_BASE_H

#include "PscFieldBase.h"

// ======================================================================
// PscInterpolatorBase

template<class G>
struct PscInterpolatorBase : PscFieldBase<interpolator_t, G>
{
  typedef PscFieldBase<interpolator_t, G> Base;
  using typename Base::Grid;
  using typename Base::Element;

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

  static PscInterpolatorBase* create(Grid *grid)
  {
    return new PscInterpolatorBase(grid);
  }
  
  static void destroy(PscInterpolatorBase* interpolator)
  {
    delete interpolator;
  }

  PscInterpolatorBase(Grid *grid)
    : Base(grid)
  {
    MALLOC_ALIGNED(arr_, grid->nv, 128);
    CLEAR(arr_, grid->nv);
  }

  ~PscInterpolatorBase()
  {
    FREE_ALIGNED(arr_);
  }
  
  Element operator[](int idx) const { return arr_[idx]; }
  Element& operator[](int idx)      { return arr_[idx]; }

private:
  using Base::arr_;

public:
  using Base::g;
};

#endif

