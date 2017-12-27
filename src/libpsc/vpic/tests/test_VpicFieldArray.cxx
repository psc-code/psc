
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

  Grid grid;
  MaterialList material_list;
  FieldArray* fa = FieldArray::create(&grid, material_list, 0.);

  test_FieldArray_methods(*fa);

  //FieldArray::destroy(fa);
}

int main(int argc, char **argv)
{
  test_VpicFieldArray();
}
