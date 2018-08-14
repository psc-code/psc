
#pragma once

#include "PscInterpolatorBase.h"

// ======================================================================
// MfieldsInterpolatorPsc

template<typename _Grid>
struct MfieldsInterpolatorPsc
{
  using Grid = _Grid;
  using Interpolator = PscInterpolatorBase<Grid>;
  using Element = typename Interpolator::Element;

  struct Patch
  {
    using Element = Element;
    
    Patch(Grid* vgrid)
      : ip_{new Interpolator{vgrid}}
    {}

    ~Patch()
    {
      delete ip_;
    }
    
    Element* data() { return ip_->data(); }
    Element  operator[](int idx) const { return (*ip_)[idx]; }
    Element& operator[](int idx)       { return (*ip_)[idx]; }

    Grid* grid() { return ip_->grid(); }

    Interpolator *ip_;
  };

  MfieldsInterpolatorPsc(Grid* vgrid)
    : patch_{vgrid}
  {}

  Patch& getPatch(int p) { return patch_; }
  
private:
  Patch patch_;
};

