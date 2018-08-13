
#pragma once

#include <fields3d.hxx>

template<typename Grid, typename HydroArray>
struct MfieldsHydroVpic_
{
  using real_t = float;
  using fields_t = fields3d<float, LayoutAOS>;

  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic_(const Grid_t& grid, Grid* vgrid)
    : grid_{grid}
  {
    assert(grid.n_patches() == 1);

    vhydro_ = new HydroArray{vgrid};
  }

  ~MfieldsHydroVpic_()
  {
    delete vhydro_;
  }

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }
  
  fields_t operator[](int p)
  {
    assert(p == 0);
    int ib[3], im[3];
    float* data = vhydro_->getData(ib, im);
    return fields_t{ib, im, N_COMP, data};
  }

  HydroArray& vhydro() { return *vhydro_; }

private:
  HydroArray* vhydro_;
  const Grid_t& grid_;
};

