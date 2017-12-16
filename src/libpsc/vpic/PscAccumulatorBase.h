
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

  using Base::Base;
};
  
// ======================================================================
// PscAccumulatorBase

template<class G>
struct PscAccumulatorBase
{
  typedef accumulator_t Element;
  typedef G Grid;
  typedef PscAccumulatorBlock<Grid> Block;

  static PscAccumulatorBase* create(Grid *grid)
  {
    return new PscAccumulatorBase(grid);
  }

  PscAccumulatorBase(Grid *grid)
    : g_(grid)
  {
    n_pipeline_ = aa_n_pipeline();
    stride_ = POW2_CEIL(g_->nv,2);
    arr_ = new Element[(n_pipeline_ + 1) * stride_]();
  }

  ~PscAccumulatorBase()
  {
    delete[] arr_;
  }

  static int aa_n_pipeline(void)
  {
    int n = serial.n_pipeline;
    if (n < thread.n_pipeline) {
      n = thread.n_pipeline;
    }
    return n;
  }

  Element* data() { return arr_; }
  
  // FIXME, not a great interface with c just another index
  Element& operator()(int c, int idx)
  {
    return arr_[c * stride_ + idx];
  }

  // FIXME, not a great interface with c just another index
  Element& operator()(int c, int i, int j, int k)
  {
    return arr_[c * stride_ + VOXEL(i,j,k, g_->nx,g_->ny,g_->nz)];
  }

  Block operator[](int c)
  {
    return Block(grid(), arr_ + c * stride_);
  }

  Grid* grid() { return g_; }

private:
  Element* arr_;

protected:
  int n_pipeline_; // Number of pipelines supported by this accumulator
  int stride_;     // Stride be each pipeline's accumulator array
  Grid* g_;
};

#endif

