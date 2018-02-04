
#include "vec3.hxx"
#include "gtest/gtest.h"

namespace {

  // ----------------------------------------------------------------------
  // Vec3Test
  
  template <typename T>
  class Vec3Test : public ::testing::Test
  {
  };

  using Vec3Types = ::testing::Types<int, float, double>;

  TYPED_TEST_CASE(Vec3Test, Vec3Types);

  TYPED_TEST(Vec3Test, ConstructorInitList)
  {
    using V3 = Vec3<TypeParam>;
    V3 v = { 1, 2, 3 };
    
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[1], 2);
    EXPECT_EQ(v[2], 3);
  }

  TYPED_TEST(Vec3Test, OperatorEqual)
  {
    using V3 = Vec3<TypeParam>;
    EXPECT_EQ((V3{ 1, 2, 3 }), (V3{ 1, 2, 3 }));
    EXPECT_NE((V3{ 1, 2, 3 }), (V3{ 1, 3, 2 }));
  }
  
  TYPED_TEST(Vec3Test, ConstructorConstValuePtr)
  {
    using V3 = Vec3<TypeParam>;
    const TypeParam arr[3] = { 1, 2, 3 };
    V3 v(arr);
    
    EXPECT_EQ(v, (V3{ 1, 2, 3 }));
  }
  
  TYPED_TEST(Vec3Test, ConstructorCopy)
  {
    using V3 = Vec3<TypeParam>;
    V3 v = { 1, 2, 3 };
    V3 v2(v);
    
    EXPECT_EQ(v2, v);
  }

  TYPED_TEST(Vec3Test, ConstructorConvertFromFloat3)
  {
    using T = TypeParam;
    using V3 = Vec3<TypeParam>;
    Vec3<float> f = { 1.2f, 2.5f, 3.8f };
    V3 v(f);
    
    EXPECT_EQ(v, (V3{ T(1.2f), T(2.5f), T(3.8f) }));
  }

  TYPED_TEST(Vec3Test, ConstructorFromScalar)
  {
    using T = TypeParam;
    using V3 = Vec3<TypeParam>;
    V3 v(T(5.f));
    
    EXPECT_EQ(v, (V3{ T(5.f), T(5.f), T(5.f) }));
  }

  TYPED_TEST(Vec3Test, DivideAssign)
  {
    using T = TypeParam;
    using V3 = Vec3<TypeParam>;
    V3 v0 = { T(1.5), T(2.), T(-3.) };
    V3 v = v0;
    V3 w = { T(1.),  T(4.), T(-1.) };

    v /= w;
    
    EXPECT_EQ(v, (V3{ v0[0]/w[0], v0[1]/w[1], v0[2]/w[2] }));
  }

  TYPED_TEST(Vec3Test, Divide)
  {
    using T = TypeParam;
    using V3 = Vec3<TypeParam>;
    V3 v = { T(1.5), T(2.), T(-3.) };
    V3 w = { T(1.),  T(4.), T(-1.) };

    EXPECT_EQ(v / w, (V3{ v[0]/w[0], v[1]/w[1], v[2]/w[2] }));
  }

  // ----------------------------------------------------------------------
  // Vec3MixedArithmetic

  TEST(Vec3MixedArithmetic, DoubleDividedByInt)
  {
    using Double3 = Vec3<double>;

    Double3 v = { 1.5, 2., -3. };
    Int3    i = { 1,  4, -1 };

    EXPECT_EQ(v / Double3(i), (Double3{ 1.5, 0.5, 3. }));
  }

}  // namespace

