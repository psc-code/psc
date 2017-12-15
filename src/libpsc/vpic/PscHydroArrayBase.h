
#ifndef PSC_HYDRO_ARRAY_BASE_H
#define PSC_HYDRO_ARRAY_BASE_H

// ======================================================================
// PscHydroArrayBase

template<class G>
struct PscHydroArrayBase
{
  typedef G Grid;
  typedef hydro_t Element;
  
  static PscHydroArrayBase* create(Grid *grid)
  {
    return new PscHydroArrayBase(grid);
  }
  
  static void destroy(PscHydroArrayBase* hydro)
  {
    delete hydro;
  }

  PscHydroArrayBase(Grid* grid)
  {
    MALLOC_ALIGNED(h, grid->nv, 128);
    g = grid;
    /* clear(); */ // now done in derived class create()
  }
  
  ~PscHydroArrayBase()
  {
    FREE_ALIGNED(h);
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
    return &h[0].jx;
  }

  Element  operator[](int idx) const { return h[idx]; }
  Element& operator[](int idx)       { return h[idx]; }

  Element *data() { return h; }

  Grid* getGrid() { return g; }

private:
  hydro_t* ALIGNED(128) h;
  
public:
  Grid* g;
};

#endif

