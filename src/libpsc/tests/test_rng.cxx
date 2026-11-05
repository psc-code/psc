#include <gtest/gtest.h>

#include <set>

#include "rng.hxx"

TEST(Rng, ProcessSeedVariesWithBaseAndRank)
{
  std::set<int> seeds;
  for (int base = 1; base <= 10; base++) {
    for (int rank = 0; rank < 8; rank++) {
      int seed = rng::detail::process_seed(rank, base);
      EXPECT_GE(seed, 0);
      seeds.insert(seed);
    }
  }
  EXPECT_EQ(seeds.size(), 80);
}

TEST(Rng, CopyOfSeededUniformSharesGenerator)
{
  rng::Uniform<double> reference(0, 1, 3);
  rng::Uniform<double> original(0, 1, 3);
  rng::Uniform<double> copy = original;
  EXPECT_EQ(original.get(), reference.get());
  EXPECT_EQ(copy.get(), reference.get());
}

// creates the process generator, so it must stay the last test
TEST(Rng, SetBaseSeedAfterConstructionThrows)
{
  rng::Uniform<double> uniform;
  EXPECT_THROW(rng::set_base_seed(1), std::logic_error);
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
