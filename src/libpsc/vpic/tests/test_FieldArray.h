
#ifndef TEST_FIELD_ARRAY_H
#define TEST_FIELD_ARRAY_H

template<typename FieldArray>
FieldArray* test_FieldArray_create()
{
  int gdims[3] = { 16, 32, 1 };
  double xl[3] = { -1., -2., -4. };
  double xh[3] = {  1. , 2. , 4. };
  int np[3] = { 1, 1, 1 };
  double dt = .05;
  double cvac = 1.;
  double eps0 = 1.;

  double dx[3];
  for (int d = 0; d < 3; d++) {
    dx[d] = (xh[d] - xl[d]) / gdims[d];
  }

  auto *grid = FieldArray::Grid::create();
  grid->setup(dx, dt, cvac, eps0);
  grid->partition_periodic_box(xl, xh, gdims, np);
  
  typename FieldArray::MaterialList material_list;
  material_list.append(material_list.create("vacuum",
					    1., 1., 1.,
					    1., 1., 1.,
					    0., 0., 0.,
					    0., 0., 0.));

  return FieldArray::create(grid, material_list, 0.);
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
  FieldArray *fa = test_FieldArray_create<FieldArray>();
  test_FieldArray_methods(*fa);
  test_FieldArray_destroy(fa);
}

#endif
