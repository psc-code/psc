
#pragma once

struct MfieldsHydroVpic
{
  using real_t = float;
  using fields_t = fields3d<float, LayoutAOS>;

  enum {
    N_COMP = 16,
  };

  MfieldsHydroVpic(const Grid_t& grid)
    : grid_{grid}
  {
    assert(grid.n_patches() == 1);

    Simulation* sim;
    psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);
    vhydro_ = new HydroArray{sim->vgrid_};
  }

  ~MfieldsHydroVpic()
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
