
#include "testing.h"
#include "test_InterpolatorBase.h"

#include "VpicGridBase.h"
#include "VpicInterpolatorBase.h"

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  typedef VpicGridBase Grid;
  typedef VpicInterpolatorBase<Grid> InterpolatorBase;

  test_InterpolatorBase<InterpolatorBase>();
}
