
#pragma once

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
};

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

