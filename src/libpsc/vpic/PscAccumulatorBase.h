
#ifndef PSC_ACCUMULATOR_BASE_H
#define PSC_ACCUMULATOR_BASE_H

#include "PscFieldBase.h"

// ======================================================================
// PscAccumulatorBlock

template<class G>
struct PscAccumulatorBlock : PscFieldBase<accumulator_t, G>
{
  typedef PscFieldBase<accumulator_t, G> Base;
  using typename Base::Element;
  using typename Base::Grid;

  PscAccumulatorBlock(Element *arr, Grid *grid)
    : Base(grid, arr)
  {
  }

 private:
  using Base::arr_;
};
  
// ======================================================================
// PscAccumulatorBase

template<class G>
struct PscAccumulatorBase : accumulator_array_t {
  typedef accumulator_t Element;
  typedef G Grid;
  typedef PscAccumulatorBlock<Grid> Block;

  static PscAccumulatorBase* create(Grid *grid)
  {
    return new PscAccumulatorBase(grid);
  }

  PscAccumulatorBase(Grid *grid)
  {
    n_pipeline = aa_n_pipeline();
    stride = POW2_CEIL(grid->nv,2);
    g = grid;
    MALLOC_ALIGNED(a, (n_pipeline+1) * stride, 128);
    CLEAR(a, (n_pipeline+1) * stride);
  }

  ~PscAccumulatorBase()
  {
    FREE_ALIGNED(a);
  }

  static int aa_n_pipeline(void)
  {
    int n = serial.n_pipeline;
    if (n < thread.n_pipeline) {
      n = thread.n_pipeline;
    }
    return n;
  }

  Element* data() { return a; }
  
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

#endif

