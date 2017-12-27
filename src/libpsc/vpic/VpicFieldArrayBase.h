
#ifndef VPIC_FIELD_ARRAY_BASE_H
#define VPIC_FIELD_ARRAY_BASE_H

#include "Field3D.h"

#include "field_advance/field_advance.h"
#define IN_sfa
#include "field_advance/standard/sfa_private.h"

#include <mrc_common.h>
#include <cassert>

// ======================================================================
// VpicSfaParams

struct VpicSfaParams: sfa_params_t
{
};

// ======================================================================
// VpicFieldArrayBase

template<class G, class ML>
struct VpicFieldArrayBase : field_array_t {
  typedef G Grid;
  typedef ML MaterialList;
  typedef field_t Element;
  typedef VpicSfaParams SfaParams;
  typedef material_coefficient_t MaterialCoefficient;
  
  enum {
    EX  = 0,
    EY  = 1,
    EZ  = 2,
    CBX = 4,
    CBY = 5,
    CBZ = 6,
    N_COMP = sizeof(field_t) / sizeof(float),
  };

  static VpicFieldArrayBase* create(Grid *grid, const MaterialList& material_list, float damp)
  {
    return static_cast<VpicFieldArrayBase*>(::new_standard_field_array(grid, material_list, damp));
  }

 public:
  float* getData(int* ib, int* im)
  {
    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    im[0] = g->nx + 2*B;
    im[1] = g->ny + 2*B;
    im[2] = g->nz + 2*B;
    ib[0] = -B;
    ib[1] = -B;
    ib[2] = -B;
    return &f[0].ex;
  }

  // These operators can be used to access the field directly,
  // though the performance isn't great, so one you use Field3D
  // when performance is important
  float operator()(int m, int i, int j, int k) const
  {
    float *f_ = &f[0].ex;
    return f_[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
  }
  
  float& operator()(int m, int i, int j, int k)
  {
    float *f_ = &f[0].ex;
    return f_[VOXEL(i,j,k, g->nx,g->ny,g->nz) * N_COMP + m];
  }

  Element  operator[](int idx) const { return f[idx]; }
  Element& operator[](int idx)       { return f[idx]; }

  Element* data() { return f; }
  
  Grid* grid() { return static_cast<Grid*>(g); }
  SfaParams& params() { return *static_cast<SfaParams*>(field_array_t::params); }
  
  // I'm keeping these for now, because I tink they're a nice interface,
  // but it doesn't scale well to other kinds of fields (as one can tell
  // from the macro use...)
#define MK_COMP_ACCESSOR(cbx)			\
  float cbx(int i, int j, int k) const		\
  {						\
    const int nx = g->nx, ny = g->ny;		\
    return f[VOXEL(i,j,k, nx,ny,nz)].cbx;	\
  }						\
						\
  float& cbx(int i, int j, int k)		\
  {						\
    const int nx = g->nx, ny = g->ny;		\
    return f[VOXEL(i,j,k, nx,ny,nz)].cbx;	\
  }

  MK_COMP_ACCESSOR(cbx)
  MK_COMP_ACCESSOR(cby)
  MK_COMP_ACCESSOR(cbz)
  MK_COMP_ACCESSOR(ex)
  MK_COMP_ACCESSOR(ey)
  MK_COMP_ACCESSOR(ez)
  
};


#endif

