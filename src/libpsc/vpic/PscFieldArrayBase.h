
#ifndef PSC_FIELD_ARRAY_BASE_H
#define PSC_FIELD_ARRAY_BASE_H

#include "grid.h"
#include "material.h"

#include "field_advance/field_advance.h"

#include <mrc_common.h>
#include <cassert>

// FIXME, this file relies on VpicFieldArray.h 
// though at least that way, the duplication is limited
#define IN_sfa
#include "field_advance/standard/sfa_private.h"
#include "VpicFieldArrayBase.h"

// the below are copies, though, skipping the kernels

inline void _field_array_ctor(field_array_t *fa, grid_t *g, const material_t *m_list, float damp)
{
  assert(g && m_list && damp >= 0.);
  MALLOC_ALIGNED( fa->f, g->nv, 128 );
  CLEAR( fa->f, g->nv );
  fa->g = g;
  fa->params = create_sfa_params( g, m_list, damp );
}

inline void _field_array_dtor(field_array_t *fa)
{
  destroy_sfa_params( (sfa_params_t *)fa->params );
  FREE_ALIGNED( fa->f );
}

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
    _field_array_ctor(this, grid->getGrid_t(), material_list, damp);
  }
  
  ~PscFieldArrayBase()
  {
    _field_array_dtor(this);
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

