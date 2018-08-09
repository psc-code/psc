
#pragma once

// ======================================================================
// VpicPushFieldsOps
//
// will only work with VpicFieldArray, though.. (maybe it should be specialized...)

template<typename FieldArray>
struct VpicPushFieldsOps
{
  static void advance_b(FieldArray& fa, double frac) { return fa.kernel->advance_b(&fa, frac); }
  static void advance_e(FieldArray& fa, double frac) { return fa.kernel->advance_e(&fa, frac); }
};

