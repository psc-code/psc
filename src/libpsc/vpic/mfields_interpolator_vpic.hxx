
#pragma once

#include "PscInterpolatorBase.h"
#include "VpicInterpolatorBase.h"

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

// ======================================================================
// MfieldsInterpolatorVpic

struct MfieldsInterpolatorVpic
{
  using Grid = VpicGridBase;
  using Element = interpolator_t;

  struct Patch
  {
    using Element = Element;
    
    Patch(Grid* vgrid)
      : ip_{new_interpolator_array(vgrid)}
    {}

    ~Patch()
    {
      delete_interpolator_array(ip_);
    }
    
    Element* data() { return ip_->i; }
    Element  operator[](int idx) const { return ip_->i[idx]; }
    Element& operator[](int idx)       { return ip_->i[idx]; }

    Grid* grid() { return static_cast<Grid*>(ip_->g); }

    interpolator_array_t *ip_;
  };

  MfieldsInterpolatorVpic(Grid* vgrid)
    : patch_{vgrid}
  {}

  Patch& getPatch(int p) { return patch_; }
  
private:
  Patch patch_;
};

