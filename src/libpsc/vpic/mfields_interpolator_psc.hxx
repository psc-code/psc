
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
// MfieldsInterpolatorPsc

template<typename _Grid>
struct MfieldsInterpolatorPsc
{
  using Grid = _Grid;
  using Element = PscInterpolatorT;

  struct Patch
  {
    using Element = Element;
    
    Patch(Grid* vgrid)
      : ip_{vgrid}
    {}

    Element* data() { return ip_.data(); }
    Element  operator[](int idx) const { return ip_[idx]; }
    Element& operator[](int idx)       { return ip_[idx]; }

    Grid* grid() { return ip_.grid(); }

    PscFieldBase<Element, Grid> ip_;
  };

  MfieldsInterpolatorPsc(Grid* vgrid)
    : patch_{vgrid}
  {}

  Patch& getPatch(int p) { return patch_; }
  
private:
  Patch patch_;
};

