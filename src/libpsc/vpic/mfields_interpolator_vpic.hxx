
#pragma once

#include "PscInterpolatorBase.h"
#include "VpicInterpolatorBase.h"

template<typename Interpolator>
struct MfieldsInterpolator_
{
  using Grid = typename Interpolator::Grid;
  
  MfieldsInterpolator_(Grid* vgrid)
  : ip_{new Interpolator{vgrid}}
  {}

  ~MfieldsInterpolator_()
  {
    delete ip_;
  }

  Interpolator& getPatch(int p) { return *ip_; }
  
  Interpolator& vip() { return *ip_; }
  
private:
  Interpolator *ip_;
};

// ======================================================================

template<typename Grid>
using MfieldsInterpolatorPsc = MfieldsInterpolator_<PscInterpolatorBase<Grid>>;

using MfieldsInterpolatorVpic = MfieldsInterpolator_<PscInterpolatorBase<VpicGridBase>>;



  
