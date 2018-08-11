
#ifndef VPIC_HYDRO_ARRAY_BASE_H
#define VPIC_HYDRO_ARRAY_BASE_H

// ======================================================================
// VpicHydroArrayBase

template<class G>
struct VpicHydroArrayBase : hydro_array_t
{
  typedef G Grid;
  typedef hydro_t Element;

  // FIXME: broken,needs to wrap
#if 0
  // ctor
  new_hydro_array(grid);
  clear_hydro_array(hydro);

  // dtor
  delete_hydro_array(hydro);
#endif
  
  float* getData(int* ib, int* im)
  {
    const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
    
    im[0] = g->nx + 2*B;
    im[1] = g->ny + 2*B;
    im[2] = g->nz + 2*B;
    ib[0] = -B;
    ib[1] = -B;
    ib[2] = -B;
    return &h[0].jx;
  }

  Element  operator[](int idx) const { return h[idx]; }
  Element& operator[](int idx)       { return h[idx]; }

  Element *data() { return h; }

  Grid* grid() { return static_cast<Grid*>(g); }
};


#endif

