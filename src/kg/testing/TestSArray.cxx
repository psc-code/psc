
#include <gtest/gtest.h>

#include <kg/SArray.h>
#include <kg/SArrayView.h>

using Real = double;
using Layout = kg::LayoutSOA;

TEST(SArray, Ctor)
{
  auto f = kg::SArray<Real, Layout>{{{1, 2, 3}, {2, 3, 4}}, 2};
}

TEST(SArray, BoundsEtc)
{
  auto f = kg::SArray<Real, Layout>{{{1, 2, 3}, {2, 3, 4}}, 2};

  EXPECT_EQ(f.ib(), Int3({1, 2, 3}));
  EXPECT_EQ(f.im(), Int3({2, 3, 4}));
  EXPECT_EQ(f.n_comps(), 2);
}

template <typename SA>
static void setSArray(SA& f)
{
  f(0, 1, 2, 3) = 123;
  f(0, 2, 2, 3) = 223;
  f(0, 1, 3, 3) = 133;
  f(0, 2, 3, 3) = 233;
  f(0, 1, 4, 3) = 143;
  f(0, 2, 4, 3) = 243;

  f(1, 1, 2, 3) = -123;
  f(1, 2, 2, 3) = -223;
  f(1, 1, 3, 3) = -133;
  f(1, 2, 3, 3) = -233;
  f(1, 1, 4, 3) = -143;
  f(1, 2, 4, 3) = -243;
}

TEST(SArray, AccessSOA)
{
  auto f = kg::SArray<Real, kg::LayoutSOA>{{{1, 2, 3}, {2, 3, 1}}, 2};
  setSArray(f);

  EXPECT_TRUE(std::equal(f.data(), f.data() + 12,
                         std::vector<Real>({123, 223, 133, 233, 143, 243, -123,
                                            -223, -133, -233, -143, -243})
                           .begin()));
}

TEST(SArray, AccessAOS)
{
  // FIXME, should also be made to work for SArray
  std::vector<Real> storage(12);
  auto f = kg::SArrayView<Real, kg::LayoutAOS>({{1, 2, 3}, {2, 3, 1}}, 2,
                                               storage.data());
  setSArray(f);

  EXPECT_TRUE(std::equal(f.data(), f.data() + 12,
                         std::vector<Real>({123, -123, 223, -223, 133, -133,
                                            233, -233, 143, -143, 243, -243})
                           .begin()));
}

TEST(SArray, data)
{
  auto f = kg::SArray<Real, Layout>{{{1, 2, 3}, {2, 3, 4}}, 2};

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
