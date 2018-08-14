
#pragma once

#include "VpicGridBase.h"

// ======================================================================
// MfieldsAccumulatorVpic

struct MfieldsAccumulatorVpic
{
  using Grid = VpicGridBase;
  using Element = accumulator_t;

  struct Block
  {
    using Grid = Grid;
    using Element = Element;
    
    Block(Element *arr, Grid *g)
      : arr_{arr}, g_{g}
    {}
    
    Element  operator[](int idx) const { return arr_[idx]; }
    Element& operator[](int idx)       { return arr_[idx]; }
    
  private:
    Element *arr_;
    Grid *g_;
  };
  
  MfieldsAccumulatorVpic(Grid* vgrid)
    : aa_{::new_accumulator_array(vgrid)}
  {}

  ~MfieldsAccumulatorVpic()
  {
    ::delete_accumulator_array(aa_);
  }

  Grid* grid() { return static_cast<Grid*>(aa_->g); }  
  int n_pipeline() { return aa_->n_pipeline; }
  int stride() { return aa_->stride; }
  
  Element* data() { return aa_->a; }
  
  // FIXME, not a great interface with arr just another index
  Element& operator()(int arr, int i, int j, int k)
  {
    return aa_->a[arr * aa_->stride + VOXEL(i,j,k, aa_->g->nx,aa_->g->ny,aa_->g->nz)];
  }

  Block operator[](int arr)
  {
    return Block(aa_->a + arr * aa_->stride, grid());
  }

  operator accumulator_array_t* () { return aa_; }
  
private:
  accumulator_array_t* aa_;
};

