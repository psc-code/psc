
#pragma once

// ======================================================================
// MfieldsHydroVpic

#include "../libpsc/vpic/VpicGridBase.h"
#include "../libpsc/vpic/VpicHydroArrayBase.h"

struct MfieldsHydroVpic
{
  using Grid = VpicGridBase;
  using HydroArray = VpicHydroArrayBase<Grid>;
  using real_t = float;
  using Element = typename HydroArray::Element;
  using fields_t = fields3d<float, LayoutAOS>;

  struct Patch
  {
    using Element = Element;

    Element* data() { return ha->data(); }
    
    Element  operator[](int idx) const { return ha->h[idx]; }
    Element& operator[](int idx)       { return ha->h[idx]; }

    Grid* grid() { return ha->grid(); }
    
    HydroArray* ha;
 };
    
  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic(const Grid_t& grid, Grid* vgrid)
    : grid_{grid}, vgrid_{vgrid}
  {
    assert(grid.n_patches() == 1);

    vhydro_ = static_cast<HydroArray*>(::new_hydro_array(vgrid));
    patch_.ha = vhydro_;
    data_ = vhydro_->getData(ib_, im_);

    ::clear_hydro_array(vhydro_);
  }

  ~MfieldsHydroVpic()
  {
    ::delete_hydro_array(vhydro_);
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return data_; }
  fields_t operator[](int p) { return {ib_, im_, N_COMP, data_}; }
  Patch& getPatch(int p) { return patch_; }

  Grid* vgrid() { return vgrid_; }

  operator HydroArray*() { return vhydro_; }

private:
  HydroArray* vhydro_;
  real_t* data_;
  Int3 ib_, im_;
  Grid* vgrid_;
  const Grid_t& grid_;
  Patch patch_;
};

