
#include "kg/Vec3.h"
#include "gtest/gtest.h"

// ----------------------------------------------------------------------
// Vec3Test

template <typename T>
class Vec3Test : public ::testing::Test
{};

using Vec3Types = ::testing::Types<int, float, double>;

TYPED_TEST_SUITE(Vec3Test, Vec3Types);

TYPED_TEST(Vec3Test, ConstructorNoArg)
{
  using V = Vec3<TypeParam>;
  V v;
  // uninitialized, nothing to check
}

TYPED_TEST(Vec3Test, ConstructorEmpty)
{
  using V = Vec3<TypeParam>;
  V v{};
  EXPECT_EQ(v[0], 0);
  EXPECT_EQ(v[1], 0);
  EXPECT_EQ(v[2], 0);
}

TYPED_TEST(Vec3Test, ConstructorInitList)
{
  using V = Vec3<TypeParam>;
  V v = {1, 2, 3};

  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
}

TYPED_TEST(Vec3Test, ConstructorInitList4)
{
  using V = kg::Vec<TypeParam, 4>;
  V v = {1, 2, 3, 4};

  EXPECT_EQ(v[0], 1);
  EXPECT_EQ(v[1], 2);
  EXPECT_EQ(v[2], 3);
  EXPECT_EQ(v[3], 4);
}

TYPED_TEST(Vec3Test, CopyCtorAssign)
{
  using V = Vec3<TypeParam>;

  V v = {1, 2, 3};
  EXPECT_EQ(v, V({1, 2, 3}));
  V v2 = v;
  V v3 = std::move(v);
  EXPECT_EQ(v2, V({1, 2, 3}));
  EXPECT_EQ(v3, V({1, 2, 3}));
  V v4, v5;
  v4 = v;
  v5 = std::move(v);
  EXPECT_EQ(v4, V({1, 2, 3}));
  EXPECT_EQ(v5, V({1, 2, 3}));
}

TYPED_TEST(Vec3Test, CopyCtorAssign4)
{
  using V = kg::Vec<TypeParam, 4>;

  V v = {1, 2, 3, 4};
  EXPECT_EQ(v, V({1, 2, 3, 4}));
  V v2 = v;
  V v3 = std::move(v);
  EXPECT_EQ(v2, V({1, 2, 3, 4}));
  EXPECT_EQ(v3, V({1, 2, 3, 4}));
  V v4, v5;
  v4 = v;
  v5 = std::move(v);
  EXPECT_EQ(v4, V({1, 2, 3, 4}));
  EXPECT_EQ(v5, V({1, 2, 3, 4}));
}

TYPED_TEST(Vec3Test, OperatorEqual)
{
  using V = Vec3<TypeParam>;
  EXPECT_EQ((V{1, 2, 3}), (V{1, 2, 3}));
  EXPECT_NE((V{1, 2, 3}), (V{1, 3, 2}));
}

TYPED_TEST(Vec3Test, OperatorEqual4)
{
  using V = kg::Vec<TypeParam, 4>;
  EXPECT_EQ((V{1, 2, 3, 4}), (V{1, 2, 3, 4}));
  EXPECT_NE((V{1, 2, 3, 4}), (V{1, 2, 3, 5}));
}

TYPED_TEST(Vec3Test, fromPointer)
{
  using V = Vec3<TypeParam>;
  const TypeParam arr[3] = {1, 2, 3};
  auto v = V::fromPointer(arr);

  EXPECT_EQ(v, (V{1, 2, 3}));
}

TYPED_TEST(Vec3Test, fromPointer4)
{
  using V = kg::Vec<TypeParam, 4>;
  const TypeParam arr[] = {1, 2, 3, 4};
  auto v = V::fromPointer(arr);

  EXPECT_EQ(v, (V{1, 2, 3, 4}));
}

TYPED_TEST(Vec3Test, ConstructorConvertFromFloat3)
{
  using T = TypeParam;
  using V = Vec3<TypeParam>;
  Vec3<float> f = {1.2f, 2.5f, 3.8f};
  V v(f);

  EXPECT_EQ(v, (V{T(1.2f), T(2.5f), T(3.8f)}));
}

TYPED_TEST(Vec3Test, ConstructorConvertFromFloat4)
{
  using T = TypeParam;
  using V = kg::Vec<TypeParam, 4>;
  kg::Vec<float, 4> f = {1.2f, 2.5f, 3.8f, 4.2f};
  V v(f);

  EXPECT_EQ(v, (V{T(1.2f), T(2.5f), T(3.8f), T(4.2f)}));
}

TYPED_TEST(Vec3Test, Negative)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};

  EXPECT_EQ(-v, (V3{-v[0], -v[1], -v[2]}));
}

TYPED_TEST(Vec3Test, AddAssign)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};
  V3 v0 = v;

  v += w;

  EXPECT_EQ(v, (V3{v0[0] + w[0], v0[1] + w[1], v0[2] + w[2]}));
}

TYPED_TEST(Vec3Test, SubtractAssign)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};
  V3 v0 = v;

  v -= w;

  EXPECT_EQ(v, (V3{v0[0] - w[0], v0[1] - w[1], v0[2] - w[2]}));
}

TYPED_TEST(Vec3Test, MultiplyAssign)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};
  V3 v0 = v;

  v *= w;

  EXPECT_EQ(v, (V3{v0[0] * w[0], v0[1] * w[1], v0[2] * w[2]}));
}

TYPED_TEST(Vec3Test, MultiplyAssignScalar)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 v0 = v;

  v *= 3;

  EXPECT_EQ(v, (V3{v0[0] * 3, v0[1] * 3, v0[2] * 3}));
}

TYPED_TEST(Vec3Test, DivideAssign)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};
  V3 v0 = v;

  v /= w;

  EXPECT_EQ(v, (V3{v0[0] / w[0], v0[1] / w[1], v0[2] / w[2]}));
}

TYPED_TEST(Vec3Test, Add)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};

  EXPECT_EQ(v + w, (V3{v[0] + w[0], v[1] + w[1], v[2] + w[2]}));
}

TYPED_TEST(Vec3Test, Subtract)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};

  EXPECT_EQ(v - w, (V3{v[0] - w[0], v[1] - w[1], v[2] - w[2]}));
}

TYPED_TEST(Vec3Test, Multiply)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};

  EXPECT_EQ(v * w, (V3{v[0] * w[0], v[1] * w[1], v[2] * w[2]}));
}

TYPED_TEST(Vec3Test, MultiplyScalar)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;

  V3 v = {T(1.5), T(2.), T(-3.)};

  EXPECT_EQ(T(2) * v, (V3{T(2) * T(1.5), T(2) * T(2.), T(2) * T(-3.)}));
}

TYPED_TEST(Vec3Test, Divide)
{
  using T = TypeParam;
  using V3 = Vec3<TypeParam>;
  V3 v = {T(1.5), T(2.), T(-3.)};
  V3 w = {T(1.), T(4.), T(-1.)};

  EXPECT_EQ(v / w, (V3{v[0] / w[0], v[1] / w[1], v[2] / w[2]}));
}

// ----------------------------------------------------------------------
// Vec3MixedArithmetic

TEST(Vec3MixedArithmetic, DoubleDividedByInt)
{
  using Double3 = Vec3<double>;

  Double3 v = {1.5, 2., -3.};
  Int3 i = {1, 4, -1};

  EXPECT_EQ(v / Double3(i), (Double3{1.5, 0.5, 3.}));
}
