
#pragma once

#include "../libpsc/vpic/PscHydroArrayBase.h"
#include <fields3d.hxx>

// ======================================================================
// MfieldsHydroPsc

template<typename Grid>
struct MfieldsHydroPsc
{
  using HydroArray = PscHydroArrayBase<Grid>;
  using real_t = float;
  using Element = typename HydroArray::Element;
  using fields_t = fields3d<real_t, LayoutAOS>;
  using Patch = PscFieldBase<Element, Grid>;
  
  enum {
    N_COMP = 16,
  };

  MfieldsHydroPsc(const Grid_t& grid, Grid* vgrid)
    : grid_{grid},
      vhydro_{new HydroArray{vgrid}},
      patch_{vgrid, vhydro_->data()}
  {
    vhydro_->getData(ib_, im_);
    assert(grid.n_patches() == 1);
  }

  ~MfieldsHydroPsc()
  {
    delete vhydro_;
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return reinterpret_cast<real_t*>(patch_.data()); }

  fields_t operator[](int p) { return {ib_, im_, N_COMP, data()}; }
  Patch& getPatch(int p) { return patch_; }
  // FIXME the above two kinds of accessing a patch worth of data needs consolidation
  
  Grid* vgrid() { return patch_.grid(); }

private:
  const Grid_t& grid_;
  HydroArray* vhydro_;
  Int3 ib_, im_;
  Patch patch_;
};

