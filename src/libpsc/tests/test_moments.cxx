
#include "gtest/gtest.h"

#include "testing.hxx"

#undef USE_CUDA

using PushParticlesTestTypes = ::testing::Types<TestConfig1vbec3dSingleYZ
#ifdef USE_CUDA
					       ,TestConfig1vbec3dCudaYZ
#endif
						>;

// FIXME, obviously this name is bad...
TYPED_TEST_CASE(PushParticlesTest, PushParticlesTestTypes);

TYPED_TEST(PushParticlesTest, Moment)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  //using Moment_n = typename TypeParam::Moment_n;
  using Moment_n = ItemMomentLoopPatches<Moment_n_1st<Mparticles, Mfields>>;
  using Bnd = typename TypeParam::Bnd;
  using particle_t = typename Mparticles::particle_t;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init particles
  particle_t prt0;
  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  auto n_prts_by_patch = std::vector<uint>{1};
  Mparticles mprts{grid};
  SetupParticles<Mparticles>::setup_particles(mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
      return prt0;
    });

  Moment_n moment_n{ppsc->grid(), ppsc->obj.comm};
  moment_n.run(mprts);
  auto& mres = moment_n.result();
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
