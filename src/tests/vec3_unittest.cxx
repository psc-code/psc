
#include "vec3.hxx"
#include "gtest/gtest.h"

namespace {

// Tests Int3 == Vec3<int>

TEST(Int3Test, ConstructorInitList)
{
  Int3 i = { 1, 2, 3 };

  EXPECT_EQ(i[0], 1);
  EXPECT_EQ(i[1], 2);
  EXPECT_EQ(i[2], 3);
}

TEST(Int3Test, OperatorEqual)
{
  EXPECT_EQ((Int3{ 1, 2, 3 }), (Int3{ 1, 2, 3 }));
  EXPECT_NE((Int3{ 1, 2, 3}), (Int3{ 1, 3, 2 }));
}

TEST(Int3Test, ConstructorConstIntPtr)
{
  const int iarr[3] = { 1, 2, 3 };
  Int3 i(iarr);

  EXPECT_EQ(i, (Int3{ 1, 2, 3 }));
}

TEST(Int3Test, ConstructorCopy)
{
  Int3 i = { 1, 2, 3 };
  Int3 j(i);

  EXPECT_EQ(j, (Int3{ 1, 2, 3 }));
}

TEST(Int3Test, ConstructorConvertFromFloat3)
{
  Vec3<float> f = { 1.2, 2.5, 3.8 };
  Int3 i(f);

  EXPECT_EQ(i, (Int3{ 1, 2, 3 }));
}


  // Tests Double3 == Vec3<double>

  using Double3 = Vec3<double>;

  TEST(Double3Test, ConstructorInitList)
  {
    Double3 v = { 1.5, 2., -3. };
    
    EXPECT_EQ(v[0], 1.5);
    EXPECT_EQ(v[1], 2.);
    EXPECT_EQ(v[2], -3.);
  }

  TEST(Double3Test, DivideAssign)
  {
    Double3 v = { 1.5, 2., -3. };
    Double3 w = { 1.,  4., -1. };

    v /= w;
    
    EXPECT_EQ(v, (Double3{ 1.5, 0.5, 3. }));
  }

  TEST(Double3Test, Divide)
  {
    Double3 v = { 1.5, 2., -3. };
    Double3 w = { 1.,  4., -1. };

    EXPECT_EQ(v / w, (Double3{ 1.5, 0.5, 3. }));
  }

}  // namespace

