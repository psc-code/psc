
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
    fa.kernel->clear_jf(&fa);
  }
  
  void clear_rhof(FieldArray& fa)
  {
    fa.kernel->clear_rhof(&fa);
  }

  void synchronize_jf(FieldArray& fa)
  {
    fa.kernel->synchronize_jf(&fa);
  }

  void compute_div_b_err(FieldArray& fa)
  {
    fa.kernel->compute_div_b_err(&fa);
  }
  
};


#endif
