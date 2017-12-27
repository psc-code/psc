
#include "testing.h"

#include "test_FieldArray.h"

#include "PscGridBase.h"
#include "PscMaterial.h"
#include "PscFieldArrayBase.h"
#include "PscFieldArrayLocalOps.h"
#include "VpicFieldArrayRemoteOps.h"
#include "PscFieldArray.h"

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

  auto *grid = new typename FieldArray::Grid;
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

void test_PscFieldArray()
{
  typedef PscGridBase Grid;
  typedef PscMaterialList MaterialList;
  typedef PscFieldArrayBase<Grid, MaterialList> FieldArrayBase;
  typedef PscFieldArrayLocalOps<FieldArrayBase> FieldArrayLocalOps;
  typedef PscFieldArrayRemoteOps<FieldArrayBase> FieldArrayRemoteOps;
  typedef PscFieldArray<FieldArrayBase, FieldArrayLocalOps, FieldArrayRemoteOps> FieldArray;

  FieldArray* fa = test_FieldArray_create<FieldArray>();

  test_FieldArray_methods(*fa);

  //FieldArray::destroy(fa);
}

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);
  test_PscFieldArray();
}
