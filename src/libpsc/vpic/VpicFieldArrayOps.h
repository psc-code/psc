
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

  void compute_div_e_err(FieldArray& fa)
  {
    fa.kernel->compute_div_e_err(&fa);
  }
  
  double compute_rms_div_b_err(FieldArray &fa)
  {
    return fa.kernel->compute_rms_div_b_err(&fa);
  }
  
  double compute_rms_div_e_err(FieldArray &fa)
  {
    return fa.kernel->compute_rms_div_e_err(&fa);
  }

  void clean_div_b(FieldArray& fa)
  {
    fa.kernel->clean_div_b(&fa);
  }

  void clean_div_e(FieldArray& fa)
  {
    fa.kernel->clean_div_e(&fa);
  }

};


#endif
