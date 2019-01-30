
#include "testing.h"
#include "test_GridBase.h"

#include "PscGridBase.h"

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  typedef PscGridBase Grid;

  test_GridBase<Grid>();
}
