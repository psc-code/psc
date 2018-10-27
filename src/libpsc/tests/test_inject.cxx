
#include "gtest/gtest.h"

#include "testing.hxx"
#include "inject_impl.hxx"

// ======================================================================
// InjectTestTarget

struct InjectTestTarget
{
  bool is_inside(double crd[3])
  {
    return (crd[1] >= 40. && crd[1] <= 120.);
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
#ifdef USE_CUDA
						TestConfig1vbec3dCudaYZ,
#endif
						TestConfig1vbec3dSingleYZ>;

TYPED_TEST_CASE(InjectTest, InjectTestTypes);

// ======================================================================
// InjectTest/Test1

TYPED_TEST(InjectTest, Test1)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using PushParticles = typename TypeParam::PushParticles;
  using Inject = typename InjectSelector<Mparticles, InjectTestTarget,
					 typename TypeParam::dim>::Inject;
  using ItemMoment = typename Inject::ItemMoment_t;
  using real_t = typename Mparticles::real_t;

  const real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  const int inject_interval = 1;
  const int inject_tau = 10;
  auto target = InjectTestTarget{};
  Inject inject{grid, inject_interval, inject_tau, 0, target}; // FIXME, can't use "auto inject = Inject{...}", though I want to

  // let's start with no particles
  Mparticles mprts{grid};

  // should have no particles
  EXPECT_EQ(mprts[0].size(), 0);

  // density should be all zero
  ItemMoment moment_n{grid, grid.comm()};
  auto& mflds_n = moment_n.result();

  moment_n.run(mprts);
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	EXPECT_EQ(mflds_n[p](0, i,j,k), 0.);
      });
  }

  // do one injection
  inject(mprts);

  // density should be equal to n_injected inside target
  real_t fac = 1. / grid.norm.cori * 
    (inject_interval * grid.dt / inject_tau) / (1. + inject_interval * grid.dt / inject_tau);
  real_t n_injected = int(fac) * grid.norm.cori;
  // FIXME, it's worth noting that given the wrong choice of interval / tau / nicell, one may end up never injecting anything because the # particles to be injected turns out to be < 1 and we always round down
  //mprintf("fac = %g, n_injected = %g\n", fac, n_injected);
  
  moment_n.run(mprts);
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	double xx[3] = {grid.patches[p].x_cc(i), grid.patches[p].y_cc(j), grid.patches[p].z_cc(k)};
	real_t n_expected = target.is_inside(xx) ? n_injected : 0.;
	EXPECT_NEAR(mflds_n[p](0, i,j,k), n_expected, eps) << "ijk " << i << ":" << j << ":" << k;
      });
  }
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
