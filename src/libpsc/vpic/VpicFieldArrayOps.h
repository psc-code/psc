
#ifndef VPIC_FIELD_ARRAY_OPS_H
#define VPIC_FIELD_ARRAY_OPS_H

// ======================================================================
// VpicFieldArrayOps

template<class FA, class FieldArrayLocalOps>
struct VpicFieldArrayOps : FieldArrayLocalOps {
  typedef FA FieldArray;
  
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

  void compute_curl_b(FieldArray& fa)
  {
    fa.kernel->compute_curl_b(&fa);
  }

  void compute_rhob(FieldArray& fa)
  {
    fa.kernel->compute_rhob(&fa);
  }

  double synchronize_tang_e_norm_b(FieldArray& fa)
  {
    return fa.kernel->synchronize_tang_e_norm_b(&fa);
  }

};

template<class B, class FieldArrayLocalOps>
struct VpicFieldArray : B, VpicFieldArrayOps<B,FieldArrayLocalOps>
{
  typedef B Base;
  typedef VpicFieldArray<B, FieldArrayLocalOps> FieldArray;

  using Base::Base;

  using Base::kernel;

  // ----------------------------------------------------------------------
  // advance
  
  void advance_b(double frac)
  {
    kernel->advance_b(this, frac);
  }

  void advance_e(double frac)
  {
    kernel->advance_e(this, frac);
  }

  // ----------------------------------------------------------------------
  // energy_f

  void energy_f(double en[6])
  {
    kernel->energy_f(en, this);
  }

  // ----------------------------------------------------------------------
  // Accumulator interfaces

  void clear_jf()
  {
    kernel->clear_jf(this);
  }
  
  void clear_rhof()
  {
    kernel->clear_rhof(this);
  }

  void synchronize_jf()
  {
    kernel->synchronize_jf(this);
  }

  void synchronize_rho()
  {
    kernel->synchronize_rho(this);
  }

};
  
#endif
