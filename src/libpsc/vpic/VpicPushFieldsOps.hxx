
#pragma once

// ======================================================================
// VpicPushFieldsOps
//
// will only work with VpicFieldArray, though.. (maybe it should be specialized...)

template<typename MfieldsState, typename FieldArray>
struct VpicPushFieldsOps
{
  static void advance_b(MfieldsState& mflds, double frac) { FieldArray* fa = mflds; return fa->kernel->advance_b(fa, frac); }
  static void advance_e(MfieldsState& mflds, double frac) { FieldArray* fa = mflds; return fa->kernel->advance_e(fa, frac); }
};

