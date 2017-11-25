
template<class FieldArray>
void
test_FieldArrayLocal_methods(FieldArray& fa)
{
  fa.local_ghost_tang_b(fa);
  fa.local_ghost_norm_e(fa);
  fa.local_ghost_div_b(fa);
  fa.local_adjust_tang_e(fa);
  fa.local_adjust_norm_b(fa);
  fa.local_adjust_div_e(fa);
  fa.local_adjust_jf(fa);
  fa.local_adjust_rhof(fa);
  fa.local_adjust_rhob(fa);
}

template<class FieldArray>
void test_FieldArray_methods(FieldArray& fa)
{
  double en[6];
  fa.energy_f(en);

  fa.advance_b(1.);
  fa.advance_e(1.);

  fa.clear_jf();
  fa.synchronize_jf();
  fa.clear_rhof();
  fa.synchronize_rho();

  fa.compute_rhob();
  fa.compute_curl_b();

  fa.synchronize_tang_e_norm_b();

  fa.compute_div_e_err();
  fa.compute_rms_div_e_err();
  fa.clean_div_e();

  fa.compute_div_b_err();
  fa.compute_rms_div_b_err();
  fa.clean_div_b();
}

