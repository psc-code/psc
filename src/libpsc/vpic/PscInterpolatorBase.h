
#ifndef PSC_INTERPOLATOR_BASE_H
#define PSC_INTERPOLATOR_BASE_H

#include "PscFieldBase.h"

struct PscInterpolatorT
{
  float ex, dexdy, dexdz, d2exdydz;
  float ey, deydz, deydx, d2eydzdx;
  float ez, dezdx, dezdy, d2ezdxdy;
  float cbx, dcbxdx;
  float cby, dcbydy;
  float cbz, dcbzdz;
  float _pad[2];  // 16-byte align
};

// ======================================================================
// PscInterpolatorBase

template<class G>
struct PscInterpolatorBase : PscFieldBase<PscInterpolatorT, G>
{
  typedef PscFieldBase<PscInterpolatorT, G> Base;
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
};

#endif

