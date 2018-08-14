
#ifndef VPIC_FIELD_ARRAY_H
#define VPIC_FIELD_ARRAY_H

// ======================================================================
// VpicCleanDivOps

template<typename MfieldsState>
struct VpicCleanDivOps
{
  static void clear_rhof(field_array_t* fa) { fa->kernel->clear_rhof(fa); }
  static void synchronize_rho(field_array_t* fa) { fa->kernel->synchronize_rho(fa); }
  static void compute_div_e_err(field_array_t* fa) { fa->kernel->compute_div_e_err(fa); }
  static double compute_rms_div_e_err(field_array_t* fa) { return fa->kernel->compute_rms_div_e_err(fa); }
  static void clean_div_e(field_array_t* fa) { fa->kernel->clean_div_e(fa); }
  static void compute_div_b_err(field_array_t* fa) { fa->kernel->compute_div_b_err(fa); }
  static double compute_rms_div_b_err(field_array_t* fa) { return fa->kernel->compute_rms_div_b_err(fa); }
  static void clean_div_b(field_array_t* fa) { fa->kernel->clean_div_b(fa); }
  static double synchronize_tang_e_norm_b(field_array_t* fa) { return fa->kernel->synchronize_tang_e_norm_b(fa); }
};

// ======================================================================
// VpicAccumulateOps

template<typename MfieldsVpic>
struct VpicAccumulateOps
{
  static void clear_jf(field_array_t* fa) { fa->kernel->clear_jf(fa); }
  static void synchronize_jf(field_array_t* fa) { fa->kernel->synchronize_jf(fa); }
  static void compute_rhob(field_array_t* fa) { fa->kernel->compute_rhob(fa); }
  static void compute_curl_b(field_array_t* fa) { fa->kernel->compute_curl_b(fa); }
};

// ======================================================================
// VpicDiagOps

template<typename MfieldsState>
struct VpicDiagOps
{
  static void energy_f(field_array_t* fa, double en[6]) { fa->kernel->energy_f(en, fa); }
};

#endif
