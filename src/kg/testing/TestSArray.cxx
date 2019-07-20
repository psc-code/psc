
#include <gtest/gtest.h>

#include <kg/SArray.h>

using Real = double;
using Layout = kg::LayoutSOA;

TEST(SArray, Ctor)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};
}

TEST(SArray, BoundsEtc)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};

  EXPECT_EQ(f.ib(), Int3({1, 2, 3}));
  EXPECT_EQ(f.im(), Int3({2, 3, 4}));
  EXPECT_EQ(f.n_comps(), 2);
  EXPECT_EQ(f.n_cells(), 2 * 3 * 4);
  EXPECT_EQ(f.size(), 2 * 3 * 4 * 2);
}

TEST(SArray, index)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};

  EXPECT_EQ(f.index(0, 1, 2, 3), 0);
  EXPECT_EQ(f.index(0, 2, 2, 3), 1);
  EXPECT_EQ(f.index(1, 1, 2, 3), 2 * 3 * 4);
}

TEST(SArray, data)
{
  auto f = kg::SArray<Real, Layout>{{1, 2, 3}, {2, 3, 4}, 2};

  EXPECT_EQ(f.data(), &f(0, 1, 2, 3));
  const auto& fc = const_cast<const kg::SArray<Real, Layout>&>(f);
  EXPECT_EQ(fc.data(), &fc(0, 1, 2, 3));
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  return rc;
}
