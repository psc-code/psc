#include <gtest/gtest.h>

#include "test_common.hxx"

#include "injector_boundary_inflow.hxx"

TEST(InjectorBoundaryInflowTest, ParticleGeneratorMaxwellianTest)
{
  int kind_idx = 0;
  Grid_t::Kind kind{1.0, 1836.0, "ion"};
  ParticleGeneratorMaxwellian::Real w = 2.0;
  ParticleGeneratorMaxwellian::Real3 mean_u{0.0, 5.0, 15.0};
  ParticleGeneratorMaxwellian::Real3 temperature{0.0, 0.0, 0.0};
  ParticleGeneratorMaxwellian::Real3 pos{1.0, 2.0, 5.0};

  ParticleGeneratorMaxwellian gen{kind_idx, kind, w, mean_u, temperature};

  auto prt = gen.get(pos, {0.0, 0.0, 0.0});

  ASSERT_EQ(prt.kind, kind_idx);
  ASSERT_EQ(prt.w, w);
  ASSERT_EQ(prt.tag, 0);
  ASSERT_EQ(prt.u, mean_u); // zero temperature => exact velocity
  ASSERT_EQ(prt.x, pos);
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
