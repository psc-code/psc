
template<class FieldArray>
void test_FieldArray_methods(FieldArray& fa)
{
  double en[6];
  fa.energy_f(en);

  fa.advance_b(1.);
}

