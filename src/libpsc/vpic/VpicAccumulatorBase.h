
#ifndef VPIC_ACCUMULATOR_BASE_H
#define VPIC_ACCUMULATOR_BASE_H

// ======================================================================
// VpicAccumulatorBlock

template<class G>
struct VpicAccumulatorBlock {
  typedef accumulator_t Element;
  typedef G Grid;

  VpicAccumulatorBlock(Element *arr, Grid *g)
    : arr_(arr), g_(g)
  {
  }

  Element operator[](int idx) const
  {
    return arr_[idx];
  }

  Element& operator[](int idx)
  {
    return arr_[idx];
  }

  //private:
  Element *arr_;
  Grid *g_;
};
  
// ======================================================================
// VpicAccumulatorBase

template<class G>
struct VpicAccumulatorBase : accumulator_array_t {
  typedef accumulator_t Element;
  typedef G Grid;
  typedef VpicAccumulatorBlock<Grid> Block;

  static VpicAccumulatorBase* create(Grid *grid)
  {
    return reinterpret_cast<VpicAccumulatorBase*>(new_accumulator_array(grid));
  }
  
  Element* data()
  {
    return a;
  }
  
  // FIXME, not a great interface with arr just another index
  Element& operator()(int arr, int idx)
  {
    return a[stride * arr + idx];
  }

  // FIXME, not a great interface with arr just another index
  Element& operator()(int arr, int i, int j, int k)
  {
    return a[arr * stride + VOXEL(i,j,k, g->nx,g->ny,g->nz)];
  }

  Block operator[](int arr)
  {
    return Block(a + arr * stride, getGrid());
  }

  Grid* getGrid()
  {
    return static_cast<Grid*>(g);
  }
};

#if 0
// ----------------------------------------------------------------------
// copied from accumulator_array.c, converted from new/delete -> ctor

static int
aa_n_pipeline(void) {
  int                       n = serial.n_pipeline;
  if( n<thread.n_pipeline ) n = thread.n_pipeline;
  return n; /* max( {serial,thread,spu}.n_pipeline ) */
}

inline void
accumulator_array_ctor(accumulator_array_t * aa, grid_t * g ) {
  if( !g ) ERROR(( "Bad grid."));
  aa->n_pipeline = aa_n_pipeline();
  aa->stride     = POW2_CEIL(g->nv,2);
  aa->g          = g;
  MALLOC_ALIGNED( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride, 128 );
  CLEAR( aa->a, (size_t)(aa->n_pipeline+1)*(size_t)aa->stride );
}

inline void
accumulator_array_dtor( accumulator_array_t * aa ) {
  if( !aa ) return;
  FREE_ALIGNED( aa->a );
}
#endif

#endif

