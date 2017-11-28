
#ifndef PSC_FIELD_ARRAY_BASE_H
#define PSC_FIELD_ARRAY_BASE_H

#include "grid.h"
#include "material.h"

#include "field_advance/field_advance.h"

#include <mrc_common.h>
#include <cassert>

// FIXME, this file relies for now on VpicFieldArray.h 
// though at least that way, the duplication is limited
// Basically, this is actually still VpicFieldArrayBase in terms
// of data layout, just some methods are adapted

#define IN_sfa
#include "field_advance/standard/sfa_private.h"
#include "VpicFieldArrayBase.h"

// ======================================================================
// PscFieldArrayBase

template<class ML>
struct PscFieldArrayBase : field_array_t
{
  typedef ML MaterialList;
  typedef field_t Element;
  
  enum {
    EX  = 0,
    EY  = 1,
    EZ  = 2,
    CBX = 4,
    CBY = 5,
    CBZ = 6,
    N_COMP = sizeof(field_t) / sizeof(float),
  };
  
  PscFieldArrayBase(Grid* grid, MaterialList material_list, float damp)
  {
    assert(grid && !material_list.empty() && damp >= 0.);
    MALLOC_ALIGNED(this->f, grid->nv, 128);
    CLEAR(this->f, grid->nv);
    this->g = grid;
    this->params = create_sfa_params(grid, material_list, damp);
  }
  
  ~PscFieldArrayBase()
  {
    destroy_sfa_params((sfa_params_t *) this->params);
    FREE_ALIGNED(this->f);
  }

  Element* data()
  {
    return f;
  }
  
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

  Element operator[](int idx) const
  {
    return f[idx];
  }
  
  Element& operator[](int idx)
  {
    return f[idx];
  }
};



#endif

