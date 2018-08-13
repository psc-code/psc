
#pragma once

#include <fields3d.hxx>

// ======================================================================
// MfieldsHydroVpic_

template<typename Grid, typename _HydroArray>
struct MfieldsHydroVpic_
{
  using real_t = float;
  using fields_t = fields3d<float, LayoutAOS>;
  using HydroArray = _HydroArray;
  
  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic_(const Grid_t& grid, Grid* vgrid)
    : grid_{grid}, vgrid_{vgrid}
  {
    assert(grid.n_patches() == 1);

    vhydro_ = new HydroArray{vgrid};
    data_ = vhydro_->getData(ib_, im_);
  }

  ~MfieldsHydroVpic_()
  {
    delete vhydro_;
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return data_; }
  fields_t operator[](int p) { return {ib_, im_, N_COMP, data_}; }

  Grid* vgrid() { return vgrid_; }

  HydroArray& vhydro() { return *vhydro_; }

private:
  HydroArray* vhydro_;
  real_t* data_;
  Int3 ib_, im_;
  Grid* vgrid_;
  const Grid_t& grid_;
};

// ======================================================================
// MfieldsHydroVpic_
//
// specialized for HydroArray == VpicHydroArray

#include "../libpsc/vpic/VpicGridBase.h"
#include "../libpsc/vpic/VpicHydroArrayBase.h"

template<>
struct MfieldsHydroVpic_<VpicGridBase, VpicHydroArrayBase<VpicGridBase>>
{
  using Grid = VpicGridBase;
  using HydroArray = VpicHydroArrayBase<Grid>;
  using real_t = float;
  using fields_t = fields3d<float, LayoutAOS>;
  
  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic_(const Grid_t& grid, Grid* vgrid)
    : grid_{grid}, vgrid_{vgrid}
  {
    assert(grid.n_patches() == 1);

    vhydro_ = static_cast<HydroArray*>(::new_hydro_array(vgrid));
    data_ = vhydro_->getData(ib_, im_);

    ::clear_hydro_array(vhydro_);
  }

  ~MfieldsHydroVpic_()
  {
    ::delete_hydro_array(vhydro_);
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return data_; }
  fields_t operator[](int p) { return {ib_, im_, N_COMP, data_}; }

  Grid* vgrid() { return vgrid_; }

  HydroArray& vhydro() { return *vhydro_; }

private:
  HydroArray* vhydro_;
  real_t* data_;
  Int3 ib_, im_;
  Grid* vgrid_;
  const Grid_t& grid_;
};

