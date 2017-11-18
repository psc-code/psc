
#ifndef VPIC_MFIELDS_H
#define VPIC_MFIELDS_H

#include "vpic_iface.h"

#include <vpic.h>

// ======================================================================
// vpic_mfields

struct vpic_mfields {
  field_array_t *field_array;

  vpic_mfields() { }

  void clear_jf();
  void synchronize_jf();
  void compute_div_b_err();
  void compute_div_e_err();
  double compute_rms_div_b_err();
  double compute_rms_div_e_err();
  void clean_div_b();
  void clean_div_e();
  void compute_curl_b();
  void clear_rhof();
  void accumulate_rho_p(vpic_mparticles *vmprts);
  void synchronize_rho();
  void compute_rhob();
  double synchronize_tang_e_norm_b();
};

struct vpic_mfields_hydro {
  hydro_array_t *hydro_array;

  vpic_mfields_hydro() { }
};

#endif

