
#pragma once

// ======================================================================
// VpicPushFieldsOps
//
// will only work with VpicFieldArray, though.. (maybe it should be specialized...)

template<typename _MfieldsState>
struct VpicPushFieldsOps
{
  using MfieldsState = _MfieldsState;
  
  static void advance_b(field_array_t* fa, double frac) { return fa->kernel->advance_b(fa, frac); }
  static void advance_e(field_array_t* fa, double frac) { return fa->kernel->advance_e(fa, frac); }
};

