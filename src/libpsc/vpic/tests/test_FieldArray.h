
#ifndef TEST_FIELD_ARRAY_H
#define TEST_FIELD_ARRAY_H

#include "test_GridBase.h"

template<typename FieldArray>
FieldArray* test_FieldArray_create(typename FieldArray::Grid *g)
{
  typename FieldArray::MaterialList material_list;
  material_list.append(material_list.create("vacuum",
					    1., 1., 1.,
					    1., 1., 1.,
					    0., 0., 0.,
					    0., 0., 0.));

  return FieldArray::create(g, material_list, 0.);
}

template<typename FieldArray>
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

template<typename FieldArray>
void test_FieldArray_destroy(FieldArray* fa)
{
  // FIXME
}

template<typename FieldArray>
void test_FieldArray()
{
  auto *g = test_GridBase_create<typename FieldArray::Grid>();
  
  FieldArray *fa = test_FieldArray_create<FieldArray>(g);
  test_FieldArray_methods(*fa);
  test_FieldArray_destroy(fa);

  test_GridBase_destroy(g);
}

#endif
