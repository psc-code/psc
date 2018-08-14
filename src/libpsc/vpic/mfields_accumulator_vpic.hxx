
#pragma once

#include "VpicGridBase.h"
#include "VpicAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorVpic

struct MfieldsAccumulatorVpic : VpicAccumulatorBase<VpicGridBase>
{
};

#include "PscAccumulatorBase.h"

// ======================================================================
// MfieldsAccumulatorPsc

template<typename Grid>
struct MfieldsAccumulatorPsc : PscAccumulatorBase<Grid>
{
};

