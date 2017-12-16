
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
  
  using Base::Base;

  static PscHydroArrayBase* create(Grid *grid)
  {
    return new PscHydroArrayBase(grid);
  }
  
  static void destroy(PscHydroArrayBase* hydro)
  {
    delete hydro;
  }

  float* getData(int* ib, int* im)
  {
    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    
    im[0] = g_->nx + 2*B;
    im[1] = g_->ny + 2*B;
    im[2] = g_->nz + 2*B;
    ib[0] = -B;
    ib[1] = -B;
    ib[2] = -B;
    return &data()[0].jx;
  }

  using Base::data;
  
protected:
  using Base::g_;
};

#endif

