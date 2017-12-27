
#include "testing.h"
#include "test_GridBase.h"

#include "PscGridBase.h"

void test_PscGridBase()
{
  typedef PscGridBase Grid;

  test_GridBase<Grid>();
}

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  test_PscGridBase();
}
