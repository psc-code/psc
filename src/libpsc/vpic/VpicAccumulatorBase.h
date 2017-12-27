
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
    return a[stride() * arr + idx];
  }

  // FIXME, not a great interface with arr just another index
  Element& operator()(int arr, int i, int j, int k)
  {
    return a[arr * stride() + VOXEL(i,j,k, g->nx,g->ny,g->nz)];
  }

  Block operator[](int arr)
  {
    return Block(a + arr * stride(), grid());
  }

  int n_pipeline() { return accumulator_array_t::n_pipeline; }
  int stride() { return accumulator_array_t::stride; }
  Grid* grid() { return static_cast<Grid*>(g); }  
};

#endif

