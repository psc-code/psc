
#ifndef VPIC_FIELD_ARRAY_LOCAL_OPS_H
#define VPIC_FIELD_ARRAY_LOCAL_OPS_H

#define IN_sfa
#include "field_advance/standard/sfa_private.h"

template<class FieldArray>
struct VpicFieldArrayLocalOps
{
  static void local_ghost_tang_b(FieldArray& fa)  { ::local_ghost_tang_b(fa.f, fa.g); }
  static void local_ghost_norm_e(FieldArray& fa)  { ::local_ghost_norm_e(fa.f, fa.g); }
  static void local_ghost_div_b(FieldArray& fa)   { ::local_ghost_div_b(fa.f, fa.g); }
  static void local_adjust_tang_e(FieldArray& fa) { ::local_adjust_tang_e(fa.f, fa.g); }
  static void local_adjust_norm_b(FieldArray& fa) { ::local_adjust_norm_b(fa.f, fa.g); }
  static void local_adjust_div_e(FieldArray& fa)  { ::local_adjust_div_e(fa.f, fa.g);  }
  static void local_adjust_jf(FieldArray& fa)     { ::local_adjust_jf(fa.f, fa.g); }
  static void local_adjust_rhof(FieldArray& fa)   { ::local_adjust_rhof(fa.f, fa.g); }
  static void local_adjust_rhob(FieldArray& fa)   { ::local_adjust_rhob(fa.f, fa.g); }
};


#endif

