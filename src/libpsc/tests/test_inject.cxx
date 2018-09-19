
#include "gtest/gtest.h"

#include "testing.hxx"
#include "../libpsc/psc_inject/psc_inject_impl.hxx"

// ======================================================================
// InjectTestTarget

struct InjectTestTarget
{
  bool is_inside(double crd[3])
  {
    return (crd[1] >= .25 && crd[1] <= .75);
  }

  void init_npt(int pop, double crd[3], struct psc_particle_npt *npt)
  {
    if (!is_inside(crd)) {
      npt->n = 0;
      return;
    }
    
    npt->n = 1.;
  }
};

// ======================================================================
// InjectTest

template<typename T>
struct InjectTest : ::testing::Test
{
  using dim = typename T::dim;
  using Mparticles = typename T::Mparticles;
  using MfieldsState = typename T::MfieldsState;
 
  const double L = 160;

  Int3 ibn = { 2, 2, 2 };

  void make_psc(const Grid_t::Kinds& kinds)
  {
    Int3 gdims = {16, 16, 16};
    if (dim::InvarX::value) { gdims[0] = 1; ibn[0] = 0; }
    if (dim::InvarY::value) { gdims[1] = 1; ibn[1] = 0; }
    if (dim::InvarZ::value) { gdims[2] = 1; ibn[2] = 0; }

    auto grid_domain = Grid_t::Domain{gdims, {L, L, L}};
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};
    
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 200;
    auto coeff = Grid_t::Normalization{norm_params};

    grid_ = Grid_t::psc_make_grid(grid_domain, grid_bc, kinds, coeff, 1., {});
  }
  
  const Grid_t& grid()
  {
    assert(grid_);
    return *grid_;
  }

  Grid_t* grid_ = {};
};

using InjectTestTypes = ::testing::Types<
#if 0//def USE_CUDA
						TestConfig1vbec3dCudaYZ,
#endif
						TestConfig1vbec3dSingleYZ>;

TYPED_TEST_CASE(InjectTest, InjectTestTypes);

// ======================================================================
// SingleParticle test

TYPED_TEST(InjectTest, SingleParticle)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  using Inject = Inject_<Mparticles, Mfields, InjectTestTarget>;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  auto inject = Inject{grid, 1, 10, 0, InjectTestTarget{}};

  // let's start with no particles
  Mparticles mprts{grid};
  mprts.inject(0, {{1., 0., 0.}, {}, 1., 0});

#if 0
  const typename Mparticles::real_t eps = 1e-5;

  // check particle
  EXPECT_EQ(mprts[0].size(), 1);
  auto prt = *mprts[0].get().begin();
  EXPECT_NEAR(prt.position()[0], 1., eps);
  EXPECT_NEAR(prt.position()[1], 0., eps);
  EXPECT_NEAR(prt.position()[2], 0., eps);
  EXPECT_NEAR(prt.w(), 1., eps);
#endif
}

// ======================================================================
// main

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
