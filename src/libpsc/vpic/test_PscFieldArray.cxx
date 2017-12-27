
#include "test_FieldArray.h"

#include "PscGridBase.h"
#include "PscMaterial.h"
#include "PscFieldArrayBase.h"
#include "PscFieldArrayLocalOps.h"
#include "VpicFieldArrayRemoteOps.h"
#include "PscFieldArray.h"

void test_PscFieldArray()
{
  typedef PscGridBase Grid;
  typedef PscMaterialList MaterialList;
  typedef PscFieldArrayBase<Grid, MaterialList> FieldArrayBase;
  typedef PscFieldArrayLocalOps<FieldArrayBase> FieldArrayLocalOps;
  typedef PscFieldArrayRemoteOps<FieldArrayBase> FieldArrayRemoteOps;
  typedef PscFieldArray<FieldArrayBase, FieldArrayLocalOps, FieldArrayRemoteOps> FieldArray;

  Grid grid;
  MaterialList material_list;
  FieldArray* fa = FieldArray::create(&grid, material_list, 0.);

  test_FieldArray_methods(*fa);

  //FieldArray::destroy(fa);
}

int main(int argc, char **argv)
{
  test_PscFieldArray();
}
