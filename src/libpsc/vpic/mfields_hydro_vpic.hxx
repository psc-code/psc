
#pragma once

// ======================================================================
// MfieldsHydroVpic

#include "../libpsc/vpic/VpicGridBase.h"

struct MfieldsHydroVpic
{
  using Grid = VpicGridBase;
  using real_t = float;
  using Element = hydro_t;
  using fields_t = fields3d<float, LayoutAOS>;

  enum {
    N_COMP = 16,
  };

  static_assert(N_COMP == sizeof(Element) / sizeof(real_t), "N_COMP doesn't match Element");

  struct Patch
  {
    using Element = Element;

    Patch(Grid* vgrid)
      : ha_{::new_hydro_array(vgrid)}
    {
      ::clear_hydro_array(ha_);
    }

    ~Patch()
    {
      ::delete_hydro_array(ha_);
    }
    
    Element* data() { return ha_->h; }
    
    Element  operator[](int idx) const { return ha_->h[idx]; }
    Element& operator[](int idx)       { return ha_->h[idx]; }

    Grid* grid() { return static_cast<Grid*>(ha_->g); }

    operator hydro_array_t*() { return ha_; }
    
  private:
    hydro_array_t* ha_;
    Grid* vgrid_;
 };
    
  MfieldsHydroVpic(const Grid_t& grid, Grid* vgrid)
    : grid_{grid},
      patch_{vgrid}
  {
    assert(grid.n_patches() == 1);

    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    im_ = { vgrid->nx + 2*B, vgrid->ny + 2*B, vgrid->nz + 2*B };
    ib_ = { -B, -B, -B };
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }

  real_t* data() { return reinterpret_cast<real_t*>(patch_.data()); }
  fields_t operator[](int p) { return {grid_, ib_, im_, N_COMP, data()}; }
  Patch& getPatch(int p) { return patch_; }

  Grid* vgrid() { return patch_.grid(); }

  operator hydro_array_t*() { return patch_; }

private:
  const Grid_t& grid_;
  Patch patch_;
  Int3 ib_, im_;
};

