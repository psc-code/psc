
#pragma once

#include "PscInterpolatorBase.h"
#include "VpicInterpolatorBase.h"

template<typename Interpolator>
struct MfieldsInterpolator_
{
  using Grid = typename Interpolator::Grid;

  struct Patch
  {
    Patch(Grid* vgrid)
      : ip_{new Interpolator{vgrid}}
    {}

    ~Patch()
    {
      delete ip_;
    }

    Interpolator *ip_;
  };

  MfieldsInterpolator_(Grid* vgrid)
    : patch_{vgrid}
  {}

  Interpolator& getPatch(int p) { return *patch_.ip_; }
  
private:
  Patch patch_;
};

// ======================================================================

template<typename Grid>
using MfieldsInterpolatorPsc = MfieldsInterpolator_<PscInterpolatorBase<Grid>>;

using MfieldsInterpolatorVpic = MfieldsInterpolator_<PscInterpolatorBase<VpicGridBase>>;



  
