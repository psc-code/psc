
#pragma once

#include "VpicGridBase.h"
#include "VpicAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorVpic

struct MfieldsAccumulatorVpic : VpicAccumulatorBase<VpicGridBase>
{
  using Base = VpicAccumulatorBase<Grid>;

  using Base::Base;
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

