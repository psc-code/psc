
#include "testing.h"

#include "test_FieldArray.h"

#include "VpicGridBase.h"
#include "VpicMaterial.h"
#include "VpicFieldArrayBase.h"
#include "VpicFieldArray.h"

void test_VpicFieldArray()
{
  typedef VpicGridBase Grid;
  typedef VpicMaterialList MaterialList;
  typedef VpicFieldArrayBase<Grid, MaterialList> FieldArrayBase;
  typedef VpicFieldArray<FieldArrayBase> FieldArray;

  FieldArray* fa = test_FieldArray_create<FieldArray>();

  test_FieldArray_methods(*fa);

  //FieldArray::destroy(fa);
}

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  test_VpicFieldArray();
}
