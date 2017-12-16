
#ifndef PSC_HYDRO_ARRAY_BASE_H
#define PSC_HYDRO_ARRAY_BASE_H

#include "PscFieldBase.h"

// ======================================================================
// PscHydroArrayBase

template<class G>
struct PscHydroArrayBase : PscFieldBase<hydro_t, G>
{
  typedef PscFieldBase<hydro_t, G> Base;
  using typename Base::Grid;
  using typename Base::Element;
  
  static PscHydroArrayBase* create(Grid *grid)
  {
    return new PscHydroArrayBase(grid);
  }
  
  static void destroy(PscHydroArrayBase* hydro)
  {
    delete hydro;
  }

  PscHydroArrayBase(Grid* grid)
    : Base(grid)
  {
    MALLOC_ALIGNED(arr_, grid->nv, 128);
    /* clear(); */ // now done in derived class create()
  }
  
  ~PscHydroArrayBase()
  {
    FREE_ALIGNED(arr_);
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
    return &arr_[0].jx;
  }

private:
  using Base::arr_;

public:
  using Base::g;
};

#endif

