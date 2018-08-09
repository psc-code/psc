
#ifndef VPIC_FIELD_ARRAY_H
#define VPIC_FIELD_ARRAY_H

// ======================================================================
// VpicCleanDivOps

template<typename FieldArray>
struct VpicCleanDivOps
{
  static void clear_rhof(FieldArray& fa) { fa.clear_rhof(); }
  static void synchronize_rho(FieldArray& fa) { fa.synchronize_rho(); }
  static void compute_div_e_err(FieldArray& fa) { fa.compute_div_e_err(); }
  static double compute_rms_div_e_err(FieldArray& fa) { return fa.compute_rms_div_e_err(); }
  static void clean_div_e(FieldArray& fa) { fa.clean_div_e(); }
  static void compute_div_b_err(FieldArray& fa) { fa.compute_div_b_err(); }
  static double compute_rms_div_b_err(FieldArray& fa) { return fa.compute_rms_div_b_err(); }
  static void clean_div_b(FieldArray& fa) { fa.clean_div_b(); }
  static double synchronize_tang_e_norm_b(FieldArray& fa) { return fa.synchronize_tang_e_norm_b(); }
};

// ======================================================================
// VpicAccumulateOps

template<typename FieldArray>
struct VpicAccumulateOps
{
  using Grid = typename FieldArray::Grid;
  
  static void clear_jf(FieldArray& fa) { fa.clear_jf(); }
  static void synchronize_jf(FieldArray& fa) { fa.synchronize_jf(); }
  static void compute_rhob(FieldArray& fa) { fa.compute_rhob(); }
  static void compute_curl_b(FieldArray& fa) { fa.compute_curl_b(); }
};

// ======================================================================
// VpicDiagOps

template<typename FieldArray>
struct VpicDiagOps
{
  static void energy_f(FieldArray& fa, double en[6]) { fa.energy_f(en); }
};

// ======================================================================
// VpicPushFieldsOps
//
// will only work with VpicFieldArray, though.. (maybe it should be specialized...)

template<typename FieldArray>
struct VpicPushFieldsOps
{
  static void advance_b(FieldArray& fa, double frac) { return fa.advance_b(frac); }
  static void advance_e(FieldArray& fa, double frac) { return fa.advance_e(frac); }
};

// ======================================================================
// VpicFieldArray

template<class B>
struct VpicFieldArray : B
{
  typedef B Base;
  typedef VpicFieldArray<B> FieldArray;
  using typename Base::Grid;
  using typename Base::MaterialList;

  using Base::kernel;

  static VpicFieldArray* create(Grid *grid, MaterialList material_list, float damp)
  {
    return static_cast<VpicFieldArray*>(Base::create(grid, material_list, damp));
  }
  
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
  // for accumulation

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

  // ----------------------------------------------------------------------
  // for initialization

  void compute_rhob()
  {
    kernel->compute_rhob(this);
  }

  void compute_curl_b()
  {
    kernel->compute_curl_b(this);
  }

  // ----------------------------------------------------------------------
  // shared face cleaning

  double synchronize_tang_e_norm_b()
  {
    return kernel->synchronize_tang_e_norm_b(this);
  }

  // ----------------------------------------------------------------------
  // div E cleaning
  
  void compute_div_e_err()
  {
    kernel->compute_div_e_err(this);
  }
  
  double compute_rms_div_e_err()
  {
    return kernel->compute_rms_div_e_err(this);
  }

  void clean_div_e()
  {
    kernel->clean_div_e(this);
  }

  // ----------------------------------------------------------------------
  // div B cleaning

  void compute_div_b_err()
  {
    kernel->compute_div_b_err(this);
  }

  double compute_rms_div_b_err()
  {
    return kernel->compute_rms_div_b_err(this);
  }
  
  void clean_div_b()
  {
    kernel->clean_div_b(this);
  }

};
  
#endif
