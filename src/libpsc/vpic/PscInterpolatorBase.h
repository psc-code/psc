
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

  using Base::Base;

  static PscInterpolatorBase* create(Grid *grid)
  {
    return new PscInterpolatorBase(grid);
  }
  
  static void destroy(PscInterpolatorBase* interpolator)
  {
    delete interpolator;
  }

public:
  using Base::g;
};

#endif

