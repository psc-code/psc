
#pragma once

#include "../libpsc/vpic/VpicGridBase.h"
#include "../libpsc/vpic/VpicMaterial.h"

#include "Field3D.h"

#include "field_advance/field_advance.h"
#define IN_sfa
#include "field_advance/standard/sfa_private.h"

#include <mrc_common.h>
#include <cassert>

// ======================================================================
// MfieldsStateVpic

struct MfieldsStateVpic
{
  using real_t = float;
  using Grid = VpicGridBase;
  using MaterialList = VpicMaterialList;
  using FieldArray = field_array_t;
  using SfaParams = sfa_params_t;

  enum {
    EX = 0, EY = 1, EZ = 2, DIV_E_ERR = 3,
    BX = 4, BY = 5, BZ = 6, DIV_B_ERR = 7,
    TCAX = 8, TCAY = 9, TCAZ = 10, RHOB = 11,
    JFX = 12, JFY = 13, JFZ = 14, RHOF = 15,
    N_COMP = 20,
  };

  struct Patch
  {
    using Element = field_t;
    
    Patch(Grid* vgrid, const MaterialList& material_list, double damp = 0.)
      : fa_{::new_standard_field_array(vgrid, material_list, damp)}
    {}

    ~Patch()
    {
      delete_field_array(fa_);
    }

    Element* data() { return fa_->f; }
    Grid* grid() { return static_cast<Grid*>(fa_->g); }
    field_t  operator[](int idx) const { return fa_->f[idx]; }
    field_t& operator[](int idx)       { return fa_->f[idx]; }

    // These operators can be used to access the field directly,
    // though the performance isn't optimal, so one should use Field3D
    // when performance is important
    static const int N_COMP = sizeof(field_t) / sizeof(float);
    
    float operator()(int m, int i, int j, int k) const
    {
      float *f = reinterpret_cast<real_t*>(fa_->f);
      return f[VOXEL(i,j,k, fa_->g->nx,fa_->g->ny,fa_->g->nz) * N_COMP + m];
    }
    
    float& operator()(int m, int i, int j, int k)
    {
      float *f = reinterpret_cast<real_t*>(fa_->f);
      return f[VOXEL(i,j,k, fa_->g->nx,fa_->g->ny,fa_->g->nz) * N_COMP + m];
    }

    SfaParams& params() { return *static_cast<SfaParams*>(fa_->params); }

    operator field_array_t* () { return fa_; }

  private:
    field_array_t* fa_;
  };
    
  using fields_t = fields3d<float, LayoutAOS>;

  MfieldsStateVpic(const Grid_t& grid, Grid* vgrid, const MaterialList& material_list, double damp = 0.)
    : grid_{grid},
      patch_{vgrid, material_list, damp}
  {
    assert(grid.n_patches() == 1);

    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    im_ = { vgrid->nx + 2*B, vgrid->ny + 2*B, vgrid->nz + 2*B };
    ib_ = { -B, -B, -B };
  }

  real_t* data() { return reinterpret_cast<real_t*>(patch_.data()); }
  fields_t operator[](int p) { return {grid(), ib_, im_, N_COMP, data()}; }
  Patch& getPatch(int p) { return patch_; }

  SfaParams& params() { return patch_.params(); }
  Grid* vgrid() { return patch_.grid(); }

  operator FieldArray*() { return patch_; }

  const Grid_t& grid() const { return grid_; }
  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return N_COMP; }
  Int3 ibn() const { return {1,1,1}; }
  
private:
  const Grid_t& grid_;
  Patch patch_;
  Int3 ib_, im_;
};

