
#include "testing.h"
#include "test_InterpolatorBase.h"

#include "PscGridBase.h"
#include "PscInterpolatorBase.h"

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  typedef PscGridBase Grid;
  typedef PscInterpolatorBase<Grid> InterpolatorBase;

  test_InterpolatorBase<InterpolatorBase>();
}
