
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

    Patch(Grid* vgrid)
      : vgrid_{vgrid}
    {}

    Element* data() { return ha->data(); }
    
    Element  operator[](int idx) const { return ha->h[idx]; }
    Element& operator[](int idx)       { return ha->h[idx]; }

    Grid* grid() { return vgrid_; }
    
    HydroArray* ha;
  private:
    Grid* vgrid_;
 };
    
  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic(const Grid_t& grid, Grid* vgrid)
    : grid_{grid},
      patch_{vgrid}
  {
    assert(grid.n_patches() == 1);

    patch_.ha = static_cast<HydroArray*>(::new_hydro_array(vgrid));
    data_ = patch_.ha->getData(ib_, im_);

    ::clear_hydro_array(patch_.ha);
  }

  ~MfieldsHydroVpic()
  {
    ::delete_hydro_array(patch_.ha);
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return data_; }
  fields_t operator[](int p) { return {ib_, im_, N_COMP, data_}; }
  Patch& getPatch(int p) { return patch_; }

  Grid* vgrid() { return patch_.grid(); }

  operator HydroArray*() { return patch_.ha; }

private:
  real_t* data_;
  Int3 ib_, im_;
  const Grid_t& grid_;
  Patch patch_;
};

