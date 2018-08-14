
#pragma once

#include "VpicGridBase.h"
#include "VpicAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorVpic

struct MfieldsAccumulatorVpic
{
  using Grid = VpicGridBase;
  using Accumulator = VpicAccumulatorBase<VpicGridBase>;
  using Element = Accumulator::Element;
  using Block = Accumulator::Block;
  
  MfieldsAccumulatorVpic(Grid* vgrid)
    : aa_{Accumulator::create(vgrid)}
  {}

  ~MfieldsAccumulatorVpic()
  {
    ::delete_accumulator_array(aa_);
  }

  Grid* grid() { return aa_->grid(); }
  int n_pipeline() const { return aa_->n_pipeline(); }
  int stride() const { return aa_->stride(); }
  
  Element* data() { return aa_->data(); }
  
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

