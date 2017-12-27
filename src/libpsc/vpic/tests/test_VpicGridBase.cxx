
#include "testing.h"
#include "test_GridBase.h"

#include "VpicGridBase.h"

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  typedef VpicGridBase Grid;

  test_GridBase<Grid>();
}
