
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

TEST(Int3Test, ConstructorConstIntPtr)
{
  const int iarr[3] = { 1, 2, 3 };
  Int3 i(iarr);

  EXPECT_EQ(i[0], 1);
  EXPECT_EQ(i[1], 2);
  EXPECT_EQ(i[2], 3);
}

TEST(Int3Test, ConstructorCopy)
{
  Int3 i = { 1, 2, 3 };
  Int3 j(i);

  EXPECT_EQ(j[0], 1);
  EXPECT_EQ(j[1], 2);
  EXPECT_EQ(j[2], 3);
}

TEST(Int3Test, ConstructorConvertType)
{
  Vec3<float> farr = { 1.2, 2.5, 3.8 };
  Int3 i(farr);

  EXPECT_EQ(i[0], 1);
  EXPECT_EQ(i[1], 2);
  EXPECT_EQ(i[2], 3);
}

}  // namespace

