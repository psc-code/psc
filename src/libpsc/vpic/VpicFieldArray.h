
#ifndef VPIC_FIELD_ARRAY_H
#define VPIC_FIELD_ARRAY_H

// ======================================================================
// VpicCleanDivOps

template<typename FieldArray>
struct VpicCleanDivOps
{
  static void clear_rhof(FieldArray& fa) { fa.kernel->clear_rhof(&fa); }
  static void synchronize_rho(FieldArray& fa) { fa.kernel->synchronize_rho(&fa); }
  static void compute_div_e_err(FieldArray& fa) { fa.kernel->compute_div_e_err(&fa); }
  static double compute_rms_div_e_err(FieldArray& fa) { return fa.kernel->compute_rms_div_e_err(&fa); }
  static void clean_div_e(FieldArray& fa) { fa.kernel->clean_div_e(&fa); }
  static void compute_div_b_err(FieldArray& fa) { fa.kernel->compute_div_b_err(&fa); }
  static double compute_rms_div_b_err(FieldArray& fa) { return fa.kernel->compute_rms_div_b_err(&fa); }
  static void clean_div_b(FieldArray& fa) { fa.kernel->clean_div_b(&fa); }
  static double synchronize_tang_e_norm_b(FieldArray& fa) { return fa.kernel->synchronize_tang_e_norm_b(&fa); }
};

// ======================================================================
// VpicAccumulateOps

template<typename FieldArray>
struct VpicAccumulateOps
{
  using Grid = typename FieldArray::Grid;
  
  static void clear_jf(FieldArray& fa) { fa.kernel->clear_jf(&fa); }
  static void synchronize_jf(FieldArray& fa) { fa.kernel->synchronize_jf(&fa); }
  static void compute_rhob(FieldArray& fa) { fa.kernel->compute_rhob(&fa); }
  static void compute_curl_b(FieldArray& fa) { fa.kernel->compute_curl_b(&fa); }
};

// ======================================================================
// VpicDiagOps

template<typename FieldArray>
struct VpicDiagOps
{
  static void energy_f(FieldArray& fa, double en[6]) { fa.kernel->energy_f(en, &fa); }
};

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

// ======================================================================
// VpicFieldArray

template<class B>
struct VpicFieldArray : B
{
  typedef B Base;
  using typename Base::Grid;
  using typename Base::MaterialList;

  static VpicFieldArray* create(Grid *grid, MaterialList material_list, float damp)
  {
    return static_cast<VpicFieldArray*>(Base::create(grid, material_list, damp));
  }
};
  
#endif
