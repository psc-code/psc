
#pragma once

#include "PscFieldBase.h"

// ======================================================================
// MfieldsAccumulatorPsc

template<typename _Grid>
struct MfieldsAccumulatorPsc
{
  using Grid = _Grid;

  // ======================================================================
  // MfieldsAccumulatorPsc::Element

  struct Element
  {
    float jx[4];   // jx0@(0,-1,-1),jx1@(0,1,-1),jx2@(0,-1,1),jx3@(0,1,1)
    float jy[4];   // jy0@(-1,0,-1),jy1@(-1,0,1),jy2@(1,0,-1),jy3@(1,0,1)
    float jz[4];   // jz0@(-1,-1,0),jz1@(1,-1,0),jz2@(-1,1,0),jz3@(1,1,0)
  };

  // ======================================================================
  // MfieldsAccumulatorPsc::Block

  struct Block : PscFieldBase<Element, Grid>
  {
    using Base = PscFieldBase<Element, Grid>;

    using Base::Base;
  };
  
  
  MfieldsAccumulatorPsc(Grid *grid)
    : g_(grid)
  {
    n_pipeline_ = aa_n_pipeline();
    stride_ = POW2_CEIL(g_->nv,2);
    arr_ = new Element[(n_pipeline_ + 1) * stride_]();
  }

  ~MfieldsAccumulatorPsc()
  {
    delete[] arr_;
  }

  static int aa_n_pipeline(void)
  {
    return 2;
#if 0
    int n = serial.n_pipeline;
    if (n < thread.n_pipeline) {
      n = thread.n_pipeline;
    }
    return n;
#endif
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
    assert(c >= 0 && c < n_pipeline_);
    return Block(grid(), arr_ + c * stride_);
  }

  int n_pipeline() { return n_pipeline_; }
  int stride() { return stride_; }
  Grid* grid() { return g_; }

private:
  Element* arr_;
  int n_pipeline_; // Number of pipelines supported by this accumulator

protected:
  int stride_;     // Stride be each pipeline's accumulator array
  Grid* g_;
};
