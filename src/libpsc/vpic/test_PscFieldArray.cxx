
#include "test_FieldArray.h"

#include "PscFieldArray.h"
#include "PscFieldArrayLocalOps.h"
#include "PscFieldArrayOps.h"

void test_PscFieldArray()
{
  typedef PscFieldArrayBase FieldArrayBase;
  typedef PscFieldArrayLocal<FieldArrayBase> FieldArrayLocal;
  typedef PscFieldArray<FieldArrayLocal> FieldArray;

  Grid grid;
  MaterialList material_list;
  FieldArray fa(&grid, material_list, 0.);

  test_FieldArray_methods(fa);
}

int main(int argc, char **argv)
{
  test_PscFieldArray();
}
