
#include <UniqueIdGenerator.h>

#include "gtest/gtest.h"

static std::unique_ptr<psc::particle::UniqueIdGenerator> uid_gen;

TEST(TestGlobalId, get)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  EXPECT_EQ((*uid_gen)(), 0 * size + rank);
  EXPECT_EQ((*uid_gen)(), 1 * size + rank);
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  uid_gen.reset(new psc::particle::UniqueIdGenerator(MPI_COMM_WORLD));
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
