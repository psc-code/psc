
#pragma once

#include "VpicGridBase.h"
#include "VpicAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorVpic

struct MfieldsAccumulatorVpic
{
  using Grid = VpicGridBase;
  using Accumulator = VpicAccumulatorBase<VpicGridBase>;
  using Element = accumulator_t;
  using Block = Accumulator::Block;
  
  MfieldsAccumulatorVpic(Grid* vgrid)
    : aa_{reinterpret_cast<Accumulator*>(::new_accumulator_array(vgrid))}
  {}

  ~MfieldsAccumulatorVpic()
  {
    ::delete_accumulator_array(aa_);
  }

  Grid* grid() { return static_cast<Grid*>(aa_->g); }  
  int n_pipeline() { return aa_->n_pipeline; }
  int stride() { return aa_->stride; }
  
  Element* data() { return aa_->a; }
  
  Element& operator()(int arr, int i, int j, int k) { return (*aa_)(arr, i, j, k); }

  Block operator[](int arr) { return (*aa_)[arr]; }

private:
  Accumulator* aa_;
};

#include "PscAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorPsc

template<typename Grid>
struct MfieldsAccumulatorPsc : PscAccumulatorBase<Grid>
{
  using Base = PscAccumulatorBase<Grid>;

  using Base::Base;
};

