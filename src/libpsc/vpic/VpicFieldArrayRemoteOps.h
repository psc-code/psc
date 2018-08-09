
#ifndef VPIC_FIELD_ARRAY_REMOTE_OPS_H
#define VPIC_FIELD_ARRAY_REMOTE_OPS_H

#define IN_sfa
#include "field_advance/standard/sfa_private.h"

template<class FieldArray>
struct VpicFieldArrayRemoteOps
{
  static void begin_remote_ghost_tang_b(FieldArray& fa) { ::begin_remote_ghost_tang_b(fa.f, fa.g); }
  static void end_remote_ghost_tang_b(FieldArray& fa)   { ::end_remote_ghost_tang_b(fa.f, fa.g);  }
  static void begin_remote_ghost_norm_e(FieldArray& fa) { ::begin_remote_ghost_norm_e(fa.f, fa.g); }
  static void end_remote_ghost_norm_e(FieldArray& fa)   { ::end_remote_ghost_norm_e(fa.f, fa.g);  }
  static void begin_remote_ghost_div_b(FieldArray& fa)  { ::begin_remote_ghost_div_b(fa.f, fa.g); }
  static void end_remote_ghost_div_b(FieldArray& fa)    { ::end_remote_ghost_div_b(fa.f, fa.g);  }
};


#endif

