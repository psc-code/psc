
#include <gtest/gtest.h>

#include <kg/SArrayView.h>

#include "fields.hxx"

using Real = double;
using Layout = kg::LayoutSOA;

template <typename SA>
static void setSArray(SA& _f)
{
  auto f = make_Fields3d<dim_xyz>(_f);
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

TEST(SArrayView, AccessSOA)
{
  std::vector<Real> storage(12);
  auto f = kg::SArrayView<Real, kg::LayoutSOA>{
    {{1, 2, 3}, {2, 3, 1}}, 2, storage.data()};
  setSArray(f);

  EXPECT_TRUE(std::equal(f.storage().data(), f.storage().data() + 12,
                         std::vector<Real>({123, 223, 133, 233, 143, 243, -123,
                                            -223, -133, -233, -143, -243})
                           .begin()));
}

TEST(SArrayView, AccessAOS)
{
  std::vector<Real> storage(12);
  auto f = kg::SArrayView<Real, kg::LayoutAOS>({{1, 2, 3}, {2, 3, 1}}, 2,
                                               storage.data());
  setSArray(f);

  EXPECT_TRUE(std::equal(f.storage().data(), f.storage().data() + 12,
                         std::vector<Real>({123, -123, 223, -223, 133, -133,
                                            233, -233, 143, -143, 243, -243})
                           .begin()));
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  return rc;
}
