
#pragma once

#include <fields3d.hxx>

// ======================================================================
// MfieldsHydroPsc

template<typename _Grid>
struct MfieldsHydroPsc
{
  using Grid = _Grid;
  using real_t = float;
  using fields_t = fields3d<real_t, LayoutAOS>;

  struct Element
  {
    float jx, jy, jz, rho; // Current and charge density => <q v_i f>, <q f>
    float px, py, pz, ke;  // Momentum and K.E. density  => <p_i f>, <m c^2 (gamma-1) f>
    float txx, tyy, tzz;   // Stress diagonal            => <p_i v_j f>, i==j
    float tyz, tzx, txy;   // Stress off-diagonal        => <p_i v_j f>, i!=j
    float _pad[2];         // 16-byte align
  };

  using Patch = PscFieldBase<Element, Grid>;
  
  enum {
    N_COMP = 16,
  };

  MfieldsHydroPsc(const Grid_t& grid, Grid* vgrid)
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
  // FIXME the above two kinds of accessing a patch worth of data needs consolidation
  
  Grid* vgrid() { return patch_.grid(); }

private:
  const Grid_t& grid_;
  Patch patch_;
  Int3 ib_, im_;
};

