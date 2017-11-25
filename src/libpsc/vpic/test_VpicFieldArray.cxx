
#include "VpicFieldArray.h"
#include "VpicFieldArrayLocalOps.h"
#include "VpicFieldArrayOps.h"

void test_VpicFieldArray()
{
  typedef VpicFieldArray<VpicFieldArrayBase, VpicFieldArrayLocalOps<VpicFieldArrayBase>> FieldArray;

  Grid grid;
  MaterialList material_list;
  FieldArray fa(&grid, material_list, 0.);
}

int main(int argc, char **argv)
{
  test_VpicFieldArray();
}
