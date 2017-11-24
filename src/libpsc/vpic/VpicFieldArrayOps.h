
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
  }

  void advance_e(FieldArray& fa, double frac)
  {
    fa.kernel->advance_e(&fa, frac);
  }

  void clear_jf(FieldArray& fa)
  {
    TIC fa.kernel->clear_jf(&fa); TOC(clear_jf, 1);
  }
  
  void clear_rhof(FieldArray& fa)
  {
    TIC fa.kernel->clear_rhof(&fa); TOC(clear_rhof, 1);
  }

};


#endif
