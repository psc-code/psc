
#ifndef VPIC_FIELD_ARRAY_H
#define VPIC_FIELD_ARRAY_H

// ======================================================================
// VpicCleanDivOps

template<typename MfieldsState, typename FieldArray>
struct VpicCleanDivOps
{
  static void clear_rhof(FieldArray* fa) { fa->kernel->clear_rhof(fa); }
  static void synchronize_rho(FieldArray* fa) { fa->kernel->synchronize_rho(fa); }
  static void compute_div_e_err(FieldArray* fa) { fa->kernel->compute_div_e_err(fa); }
  static double compute_rms_div_e_err(FieldArray* fa) { return fa->kernel->compute_rms_div_e_err(fa); }
  static void clean_div_e(FieldArray* fa) { fa->kernel->clean_div_e(fa); }
  static void compute_div_b_err(FieldArray* fa) { fa->kernel->compute_div_b_err(fa); }
  static double compute_rms_div_b_err(FieldArray* fa) { return fa->kernel->compute_rms_div_b_err(fa); }
  static void clean_div_b(FieldArray* fa) { fa->kernel->clean_div_b(fa); }
  static double synchronize_tang_e_norm_b(FieldArray* fa) { return fa->kernel->synchronize_tang_e_norm_b(fa); }
};

// ======================================================================
// VpicAccumulateOps

template<typename MfieldsVpic, typename FieldArray>
struct VpicAccumulateOps
{
  using Grid = typename MfieldsVpic::Grid;
  
  static void clear_jf(FieldArray* fa) { fa->kernel->clear_jf(fa); }
  static void synchronize_jf(FieldArray* fa) { fa->kernel->synchronize_jf(fa); }
  static void compute_rhob(FieldArray* fa) { fa->kernel->compute_rhob(fa); }
  static void compute_curl_b(FieldArray* fa) { fa->kernel->compute_curl_b(fa); }
};

// ======================================================================
// VpicDiagOps

template<typename MfieldsState, typename FieldArray>
struct VpicDiagOps
{
  static void energy_f(FieldArray* fa, double en[6]) { fa->kernel->energy_f(en, fa); }
};

#endif
