
#ifndef VPIC_FIELD_ARRAY_OPS_H
#define VPIC_FIELD_ARRAY_OPS_H

// ======================================================================
// VpicFieldArrayOps

template<class FA>
struct VpicFieldArrayOps {
  typedef FA FieldArray;
  
  void advance_b(FieldArray& fa, double frac)
  {
    fa.kernel->advance_b(&fa, frac);
    //advanceB(frac);
  }
};


#endif
