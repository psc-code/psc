
#include <kg/io.h>

#include <gtest/gtest.h>

TEST(IO, WriteRead)
{
  auto io = kg::io::IOAdios2{};

  {
    auto writer = io.open("test.bp", kg::io::Mode::Write);
    writer.beginStep(kg::io::StepMode::Append);
    writer.putLocal("var_double", 99.);
    writer.endStep();
    writer.close();
  }

  {
    auto reader = io.open("test.bp", kg::io::Mode::Read);
    reader.beginStep(kg::io::StepMode::Read);
    double dbl;
    reader.getLocal("var_double", dbl);
    reader.endStep();
    reader.close();
    EXPECT_EQ(dbl, 99.);
  }
}

TEST(IO, WriteReadAttr)
{
  auto io = kg::io::IOAdios2{};

  {
    auto writer = io.open("test.bp", kg::io::Mode::Write);
    writer.beginStep(kg::io::StepMode::Append);
    writer.put("attr_double", 99.);
    writer.endStep();
    writer.close();
  }

  {
    auto reader = io.open("test.bp", kg::io::Mode::Read);
    reader.beginStep(kg::io::StepMode::Read);
    double dbl;
    reader.get("attr_double", dbl);
    reader.endStep();
    reader.close();
    EXPECT_EQ(dbl, 99.);
  }
}

struct Custom
{
  int i;
  double d;
};

template <>
class kg::io::Descr<Custom>
{
public:
  static void put(kg::io::Engine& writer, const Custom& c,
                  const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    writer.put("i", c.i, launch);
    writer.putLocal("d", c.d, launch);
  }

  static void get(kg::io::Engine& reader, Custom& c,
                  const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    reader.get("i", c.i, launch);
    reader.getLocal("d", c.d, launch);
  }
};

TEST(IO, WriteReadCustom)
{
  auto io = kg::io::IOAdios2{};

  {
    auto writer = io.open("test.bp", kg::io::Mode::Write);
    writer.beginStep(kg::io::StepMode::Append);
    auto c = Custom{3, 99.};
    writer.put("var_custom", c);
    writer.endStep();
    writer.close();
  }

  {
    auto reader = io.open("test.bp", kg::io::Mode::Read);
    reader.beginStep(kg::io::StepMode::Read);
    auto c = Custom{};
    reader.get("var_custom", c);
    reader.endStep();
    reader.close();
    EXPECT_EQ(c.i, 3);
    EXPECT_EQ(c.d, 99.);
  }
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
