
#include "range.hxx"
#include "gtest/gtest.h"

namespace {

  // ----------------------------------------------------------------------
  // RangeTest
  
  template <typename T>
  class RangeTest : public ::testing::Test
  {
  };

  using RangeTypes = ::testing::Types<int>;

  TYPED_TEST_CASE(RangeTest, RangeTypes);

  TYPED_TEST(RangeTest, OneParam)
  {
    using T = TypeParam;

    std::vector<T> vals;
    for (auto i : range(T{5})) {
      vals.push_back(i);
    }
    EXPECT_EQ(vals, std::vector<T>({ 0, 1, 2, 3, 4 }));
  }

  TYPED_TEST(RangeTest, TwoParams)
  {
    using T = TypeParam;

    std::vector<T> vals;
    for (auto i : range(T{10}, T{15})) {
      vals.push_back(i);
    }
    EXPECT_EQ(vals, std::vector<T>({ 10, 11, 12, 13, 14 }));
  }

  TYPED_TEST(RangeTest, ThreeParams)
  {
    using T = TypeParam;

    std::vector<T> vals;
    for (auto i : range(T{10}, T{15}, T{2})) {
      vals.push_back(i);
    }
    EXPECT_EQ(vals, std::vector<T>({ 10, 12, 14 }));

    vals.clear();
    for (auto i : range(T{10}, T{14}, T{2})) {
      vals.push_back(i);
    }
    EXPECT_EQ(vals, std::vector<T>({ 10, 12 }));
}

}  // namespace

