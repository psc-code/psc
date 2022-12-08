
#include <gtest/gtest.h>

#include <psc/atomic.hxx>

TEST(Atomic, Convert)
{
  psc::atomic<float> x = 2.;
  EXPECT_EQ(x, 2.);

  x += 3.;
  EXPECT_EQ(x, 5.);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  return rc;
}
