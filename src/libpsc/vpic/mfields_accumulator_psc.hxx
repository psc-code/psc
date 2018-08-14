
#pragma once

#include "PscAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorPsc

template<typename Grid>
struct MfieldsAccumulatorPsc : PscAccumulatorBase<Grid>
{
  using Base = PscAccumulatorBase<Grid>;

  using Base::Base;
};

