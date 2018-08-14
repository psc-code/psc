
#pragma once

// ======================================================================
// MfieldsHydroVpic

#include "../libpsc/vpic/VpicGridBase.h"
#include "../libpsc/vpic/VpicHydroArrayBase.h"

struct MfieldsHydroVpic
{
  using Grid = VpicGridBase;
  using real_t = float;
  using Element = hydro_t;
  using fields_t = fields3d<float, LayoutAOS>;

  struct Patch
  {
    using Element = Element;

    Patch(Grid* vgrid)
      : vgrid_{vgrid}
    {}

    Element* data() { return ha->h; }
    
    Element  operator[](int idx) const { return ha->h[idx]; }
    Element& operator[](int idx)       { return ha->h[idx]; }

    Grid* grid() { return vgrid_; }
    
    hydro_array_t* ha;
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

    patch_.ha = ::new_hydro_array(vgrid);

    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    im_ = { vgrid->nx + 2*B, vgrid->ny + 2*B, vgrid->nz + 2*B };
    ib_ = { -B, -B, -B };
    data_ = &patch_.ha->h[0].jx;

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

  operator hydro_array_t*() { return patch_.ha; }

private:
  real_t* data_;
  Int3 ib_, im_;
  const Grid_t& grid_;
  Patch patch_;
};

