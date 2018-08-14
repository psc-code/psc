
#pragma once

#include "../libpsc/vpic/VpicGridBase.h"
#include "../libpsc/vpic/VpicMaterial.h"
#include "../libpsc/vpic/VpicFieldArrayBase.h"

// ======================================================================
// MfieldsStateVpic

struct MfieldsStateVpic
{
  using real_t = float;
  using Grid = VpicGridBase;
  using MaterialList = VpicMaterialList;
  using FieldArray = VpicFieldArrayBase<Grid, MaterialList>;

  enum {
    EX = 0, EY = 1, EZ = 2, DIV_E_ERR = 3,
    BX = 4, BY = 5, BZ = 6, DIV_B_ERR = 7,
    TCAX = 8, TCAY = 9, TCAZ = 10, RHOB = 11,
    JFX = 12, JFY = 13, JFZ = 14, RHOF = 15,
    N_COMP = 20,
  };

  using fields_t = fields3d<float, LayoutAOS>;
  using Patch = FieldArray;

  MfieldsStateVpic(const Grid_t& grid, Grid* vgrid, const MaterialList& material_list, double damp = 0.)
    : grid_{grid}
  {
    assert(grid.n_patches() == 1);

    vmflds_fields_ = FieldArray::create(vgrid, material_list, damp);
  }

  const Grid_t& grid() const { return grid_; }
  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }
  Int3 ibn() const { return {1,1,1}; }
  
  fields_t operator[](int p)
  {
    assert(p == 0);
    int ib[3], im[3];
    float* data = vmflds_fields_->getData(ib, im);
    return {ib, im, N_COMP, data};
  }

  Patch& getPatch(int p) { return *vmflds_fields_; }
  Grid* vgrid() { return vmflds_fields_->grid(); }

  operator FieldArray*() { return vmflds_fields_; }

private:
  FieldArray* vmflds_fields_;
  const Grid_t& grid_;
};

