
#include <gtest/gtest.h>

#include "../vpic/PscRng.h"

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

TEST(Rng, RngPool)
{
  RngPool rngpool;
  Rng *rng = rngpool[0];

  for (int i = 0; i < 100; ++i) {
    double r = rng->uniform(0., 1000.);
    EXPECT_GE(r, 0.);
    EXPECT_LE(r, 1000.);
  }
}

