
#include <kg/io/IOAdios2.h>

#include <gtest/gtest.h>

TEST(IOAdios2, CtorDtor) { auto io = kg::io::IOAdios2{}; }

TEST(IOAdios2, OpenWrite)
{
  auto io = kg::io::IOAdios2{};
  auto file = io.openFile("test1", kg::io::Mode::Write);
}

TEST(IOAdios2, OpenReadMissingFile)
{
  auto io = kg::io::IOAdios2{};
  EXPECT_THROW(io.openFile("test_missing", kg::io::Mode::Read),
               std::ios_base::failure);
}

TEST(IOAdios2, OpenWriteThenRead)
{
  auto io = kg::io::IOAdios2{};
  {
    auto file = io.openFile("test2.bp", kg::io::Mode::Write);
  }
  {
    auto file = io.openFile("test2.bp", kg::io::Mode::Read);
  }
}

TEST(IOAdios2, FilePutGetVariable)
{
  auto io = kg::io::IOAdios2{};
  {
    auto file = io.openFile("test3.bp", kg::io::Mode::Write);
    file.beginStep(kg::io::StepMode::Append);
    auto dbl = std::vector<double>{1., 2., 3., 4., 5.};
    file.putVariable("dbl", dbl.data(), kg::io::Mode::NonBlocking, {5},
                     {{0}, {5}}, {});
    file.putVariable("dbl2", dbl.data(), kg::io::Mode::NonBlocking, {5},
                     {{0}, {2}}, {});
    file.putVariable("dbl2", dbl.data() + 3, kg::io::Mode::NonBlocking, {5},
                     {{3}, {2}}, {});
    file.performPuts();
    file.endStep();
  }
  {
    auto file = io.openFile("test3.bp", kg::io::Mode::Read);
    file.beginStep(kg::io::StepMode::Read);

    auto shape = file.shapeVariable("dbl");
    EXPECT_EQ(shape, kg::io::Dims{5});
    auto dbl = std::vector<double>(shape[0]);
    file.getVariable("dbl", dbl.data(), kg::io::Mode::NonBlocking, {{0}, {5}},
                     {});

    auto shape2 = file.shapeVariable("dbl2");
    EXPECT_EQ(shape2, kg::io::Dims{5});
    auto dbl2 = std::vector<double>(shape[0]);
    file.getVariable("dbl2", dbl2.data(), kg::io::Mode::NonBlocking, {{0}, {5}},
                     {});
    file.performGets();
    file.endStep();

    EXPECT_EQ(dbl, (std::vector<double>{1., 2., 3., 4., 5.}));
    EXPECT_EQ(dbl2, (std::vector<double>{1., 2., 0., 4., 5.}));
  }
}

TEST(IOAdios2, FilePutGetAttribute)
{
  auto io = kg::io::IOAdios2{};
  {
    auto file = io.openFile("test4.bp", kg::io::Mode::Write);
    file.beginStep(kg::io::StepMode::Append);
    auto dbl = std::vector<double>{1., 2., 3., 4., 5.};
    file.putAttribute("attr_dbl", dbl.data(), dbl.size());
    file.endStep();
  }
  {
    auto file = io.openFile("test4.bp", kg::io::Mode::Read);
    file.beginStep(kg::io::StepMode::Read);
    auto size = file.sizeAttribute("attr_dbl");
    auto dbl = std::vector<double>(size);
    file.getAttribute("attr_dbl", dbl.data());
    EXPECT_EQ(dbl, (std::vector<double>{1., 2., 3., 4., 5.}));
    file.endStep();
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
