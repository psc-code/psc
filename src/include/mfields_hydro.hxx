
#pragma once

#include <fields3d.hxx>

// ======================================================================
// MfieldsHydroVpic_

template<typename Grid, typename _HydroArray>
struct MfieldsHydroVpic_
{
  using real_t = float;
  using HydroArray = _HydroArray;
  using Element = typename HydroArray::Element;
  using fields_t = fields3d<float, LayoutAOS>;
  using Patch = PscFieldBase<Element, Grid>;
  
  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic_(const Grid_t& grid, Grid* vgrid)
    : grid_{grid},
      vgrid_{vgrid},
      vhydro_{new HydroArray{vgrid}},
      data_{vhydro_->getData(ib_, im_)},
      patch_{vgrid, reinterpret_cast<Element*>(data_)}
  {
    assert(grid.n_patches() == 1);
  }

  ~MfieldsHydroVpic_()
  {
    delete vhydro_;
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return data_; }

  fields_t operator[](int p) { return {ib_, im_, N_COMP, data_}; }
  Patch& getPatch(int p) { return patch_; }
  // FIXME the above two kinds of accessing a patch worth of data needs consolidation
  
  Grid* vgrid() { return vgrid_; }

private:
  const Grid_t& grid_;
  Grid* vgrid_;
  HydroArray* vhydro_;
  real_t* data_;
  Int3 ib_, im_;
  Patch patch_;
};

// ======================================================================
// MfieldsHydroVpic

#include "../libpsc/vpic/VpicGridBase.h"
#include "../libpsc/vpic/VpicHydroArrayBase.h"

struct MfieldsHydroVpic
{
  using Grid = VpicGridBase;
  using HydroArray = VpicHydroArrayBase<Grid>;
  using real_t = float;
  using fields_t = fields3d<float, LayoutAOS>;
  using Patch = HydroArray;
    
  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic(const Grid_t& grid, Grid* vgrid)
    : grid_{grid}, vgrid_{vgrid}
  {
    assert(grid.n_patches() == 1);

    vhydro_ = static_cast<HydroArray*>(::new_hydro_array(vgrid));
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
  Patch& getPatch(int p) { return *vhydro_; }

  Grid* vgrid() { return vgrid_; }

  operator HydroArray*() { return vhydro_; }

private:
  HydroArray* vhydro_;
  real_t* data_;
  Int3 ib_, im_;
  Grid* vgrid_;
  const Grid_t& grid_;
};

