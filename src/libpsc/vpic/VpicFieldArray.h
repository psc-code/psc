
#ifndef VPIC_FIELD_ARRAY_H
#define VPIC_FIELD_ARRAY_H

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
