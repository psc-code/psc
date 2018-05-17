
#include "grid.hxx"
#include "fields.hxx"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "push_particles.hxx"
#include "../libpsc/psc_push_particles/push_config.hxx"
#include "../libpsc/psc_push_particles/push_dispatch.hxx"
#include "../libpsc/psc_push_particles/1vb/push_particles_1vbec_single.hxx"
#include "setup_fields.hxx"
#include "../libpsc/psc_checks/checks_impl.hxx"
#include "setup_particles.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/push_particles_cuda_impl.hxx"
#include "../libpsc/cuda/setup_fields_cuda.hxx"
#include "../libpsc/cuda/setup_particles_cuda.hxx"
#endif

#include "gtest/gtest.h"

// Rng hackiness

#include "../vpic/PscRng.h"

using Rng = PscRng;
using RngPool = PscRngPool<Rng>;

template<class Mparticles>
struct Getter
{
  Getter(Mparticles& mprts)
    : mprts_(mprts) {}

  typename Mparticles::patch_t& operator[](int p) { return mprts_[p]; }

private:
  Mparticles& mprts_;
};

template<>
struct Getter<MparticlesCuda<BS144>>
{
  using Mparticles = MparticlesCuda<BS144>;

  Getter(Mparticles& mprts)
    : mprts_(mprts.template get_as<MparticlesSingle>()) {}

  typename MparticlesSingle::patch_t& operator[](int p) { return mprts_[p]; }
  

private:
  MparticlesSingle& mprts_;
};

template<class Mparticles>
Getter<Mparticles> make_getter(Mparticles& mprts)
{
  return Getter<Mparticles>(mprts);
}

// ======================================================================
// class PushParticlesTest

template<typename T>
struct PushParticlesTest : ::testing::Test
{
  const double L = 160;

  ~PushParticlesTest()
  {
    ppsc = NULL; // FIXME, should use psc_destroy(ppsc), or really, get rid of ppsc...
  }

  void make_psc(const Grid_t::Kinds& kinds)
  {
    Int3 gdims = {16, 16, 16};
    if (T::dim::InvarX::value) gdims[0] = 1;
    if (T::dim::InvarY::value) gdims[1] = 1;
    if (T::dim::InvarZ::value) gdims[2] = 1;

    auto grid_domain = Grid_t::Domain{gdims, {L, L, L}};
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};
    
    auto psc = psc_create(MPI_COMM_WORLD); // to create ppsc, mostly
    psc_default_dimensionless(psc);
    psc_setup_coeff(psc);
    
    psc_setup_domain(psc, grid_domain, grid_bc, kinds);
    
    psc->dt = psc->grid_->dt = 1.;
  }
  
  const Grid_t& grid()
  {
    assert(ppsc);
    return ppsc->grid();
  }
};

template<typename DIM, typename PUSHP, typename ORDER>
struct TestConfig
{
  using dim = DIM;
  using PushParticles = PUSHP;
  using Mparticles = typename PushParticles::Mparticles;
  using Mfields = typename PushParticles::Mfields;
  using Checks = Checks_<Mparticles, Mfields, ORDER>;
};

using TestConfig2ndDouble = TestConfig<dim_xyz,
				       PushParticles__<Config2ndDouble<dim_xyz>>,
				       checks_order_2nd>;
using TestConfig2ndDoubleYZ = TestConfig<dim_yz,
				       PushParticles__<Config2ndDouble<dim_xyz>>,
				       checks_order_2nd>;
using TestConfig2ndSingle = TestConfig<dim_xyz,
				       PushParticles__<Config2nd<MparticlesSingle, MfieldsSingle, dim_xyz>>,
				       checks_order_2nd>;
using TestConfig1vbec3dSingle = TestConfig<dim_xyz,
					   PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsSingle, dim_xyz>>,
					   checks_order_1st>;
using TestConfig1vbec3dSingleYZ = TestConfig<dim_yz,
					     PushParticles1vb<Config1vbecSplit<MparticlesSingle, MfieldsSingle, dim_yz>>,
					     checks_order_1st>;

#ifdef USE_CUDA
using Config1vbec3d = PushParticlesConfig<BS144, opt_ip_1st_ec, DepositVb3d, CurrmemShared>;
using TestConfig1vbec3dCudaYZ = TestConfig<dim_xyz,
					   PushParticlesCuda<Config1vbec3d>,
					   checks_order_1st>;
#endif

using PushParticlesTestTypes = ::testing::Types<TestConfig2ndDoubleYZ,
						TestConfig1vbec3dSingleYZ,
#ifdef USE_CUDA
						//						TestConfig1vbec3dCudaYZ,
#endif
						TestConfig2ndDouble,
						TestConfig2ndSingle,
						TestConfig1vbec3dSingle>;

TYPED_TEST_CASE(PushParticlesTest, PushParticlesTestTypes);

// ======================================================================
// SingleParticle test

TYPED_TEST(PushParticlesTest, SingleParticle)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init particle
  auto n_prts_by_patch = std::vector<uint>{1};

  Mparticles mprts{grid};
  SetupParticles<Mparticles>::setup_particles(mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
      typename Mparticles::particle_t prt{};
      prt.xi = 1.;
      prt.yi = 0.;
      prt.zi = 0.;
      prt.qni_wni_ = 1.;
      return prt;
    });

  // check particle
  for (auto& prt : make_getter(mprts)[0]) {
    EXPECT_NEAR(prt.xi, 1., eps);
    EXPECT_NEAR(prt.yi, 0., eps);
    EXPECT_NEAR(prt.zi, 0., eps);
    EXPECT_NEAR(prt.qni_wni_, 1., eps);
  }
}

// ======================================================================
// SingleParticlePushp test

TYPED_TEST(PushParticlesTest, SingleParticlePushp)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init fields
  auto mflds = Mfields{grid, NR_FIELDS, {2, 2, 2}};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  // init particle
  auto n_prts_by_patch = std::vector<uint>{1};

  Mparticles mprts{grid};
  SetupParticles<Mparticles>::setup_particles(mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
      typename Mparticles::particle_t prt{};
      prt.xi = this->L/2; prt.yi = this->L/2; prt.zi = this->L/2;
      prt.pxi = 0.; prt.pyi = 0.; prt.pzi = 0.;
      prt.qni_wni_ = 1.;
      return prt;
    });

  PushParticles pushp_;
  pushp_.push_mprts(mprts, mflds);
  
  mprts.dump("prts.asc");

  for (auto& prt : make_getter(mprts)[0]) {
    EXPECT_NEAR(prt.pxi, 1., eps);
    EXPECT_NEAR(prt.pyi, 2., eps);
    EXPECT_NEAR(prt.pzi, 3., eps);
    EXPECT_NEAR(prt.qni_wni_, 1., eps);
  }
}

// ======================================================================
// PushParticlesTest2 is the same, but won't test cuda (because checks doesn't work)

template<typename T>
struct PushParticlesTest2 : PushParticlesTest<T>
{};

using PushParticlesTest2Types = ::testing::Types<TestConfig2ndDouble,
						 TestConfig2ndSingle,
						 TestConfig1vbec3dSingle>;

TYPED_TEST_CASE(PushParticlesTest2, PushParticlesTest2Types);

// ======================================================================
// Accel test

TYPED_TEST(PushParticlesTest2, Accel)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  using Checks = typename TypeParam::Checks;

  const int n_prts = 131;
  const int n_steps = 10;
  const typename Mparticles::real_t eps = 1e-5;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init fields
  auto mflds = Mfields{grid, NR_FIELDS, {2, 2, 2}};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case EX: return 1.;
      case EY: return 2.;
      case EZ: return 3.;
      default: return 0.;
      }
    });

  // init particles
  RngPool rngpool;
  Rng *rng = rngpool[0];
  auto n_prts_by_patch = std::vector<uint>{n_prts};

  Mparticles mprts{grid};
  SetupParticles<Mparticles>::setup_particles(mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
      typename Mparticles::particle_t prt{};
      prt.xi = rng->uniform(0, this->L);
      prt.yi = rng->uniform(0, this->L);
      prt.zi = rng->uniform(0, this->L);
      prt.qni_wni_ = 1.;
      prt.pxi = 0.;
      prt.pyi = 0.;
      prt.pzi = 0.;
      prt.kind_ = 0;
      return prt;
    });

  // run test
  PushParticles pushp_;
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    //checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    //checks_.continuity_after_particle_push(mprts, mflds);

    for (auto& prt : make_getter(mprts)[0]) {
      EXPECT_NEAR(prt.pxi, 1*(n+1), eps);
      EXPECT_NEAR(prt.pyi, 2*(n+1), eps);
      EXPECT_NEAR(prt.pzi, 3*(n+1), eps);
      prt.xi = this->L/2;
      prt.yi = this->L/2;
      prt.zi = this->L/2;
    }
  }
}

// ======================================================================
// Cyclo test

TYPED_TEST(PushParticlesTest2, Cyclo)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using PushParticles = typename TypeParam::PushParticles;
  using Checks = typename TypeParam::Checks;

  const int n_prts = 131;
  const int n_steps = 64;
  // the errors here are (substantial) truncation error, not
  // finite precision, and they add up
  // (but that's okay, if a reminder that the 6th order Boris would
  //  be good)
  const typename Mparticles::real_t eps = 1e-2;

  auto kinds = Grid_t::Kinds{Grid_t::Kind(2., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init fields
  auto mflds = Mfields{grid, NR_FIELDS, {2, 2, 2}};
  SetupFields<Mfields>::set(mflds, [&](int m, double crd[3]) {
      switch (m) {
      case HZ: return 2. * M_PI / n_steps;
      default: return 0.;
      }
    });

  // init particles
  RngPool rngpool;
  Rng *rng = rngpool[0];
  auto n_prts_by_patch = std::vector<uint>{n_prts};

  Mparticles mprts{grid};
  mprts.reserve_all(n_prts_by_patch.data());
  for (int n = 0; n < n_prts_by_patch[0]; n++) {
    typename Mparticles::particle_t prt{};
    prt.xi = rng->uniform(0, this->L);
    prt.yi = rng->uniform(0, this->L);
    prt.zi = rng->uniform(0, this->L);
    prt.pxi = 1.; // gamma = 2
    prt.pyi = 1.;
    prt.pzi = 1.;
    prt.qni_wni_ = rng->uniform(0, 1.);;
    prt.kind_ = 0;
    mprts[0].push_back(prt);
  }

  // run test
  PushParticles pushp_;
  ChecksParams checks_params{};
  checks_params.continuity_threshold = 1e-10;
  checks_params.continuity_verbose = false;
  Checks checks_{grid, MPI_COMM_WORLD, checks_params};
  for (int n = 0; n < n_steps; n++) {
    //checks_.continuity_before_particle_push(mprts);
    pushp_.push_mprts(mprts, mflds);
    //checks_.continuity_after_particle_push(mprts, mflds);

    double ux = (cos(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 cos(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uy = (sin(2*M_PI*(0.125*n_steps-(n+1))/(double)n_steps) /
		 sin(2*M_PI*(0.125*n_steps)      /(double)n_steps));
    double uz = 1.;
    for (auto& prt : make_getter(mprts)[0]) {
      EXPECT_NEAR(prt.pxi, ux, eps);
      EXPECT_NEAR(prt.pyi, uy, eps);
      EXPECT_NEAR(prt.pzi, uz, eps);
      prt.xi = this->L/2;
      prt.yi = this->L/2;
      prt.zi = this->L/2;
    }
  }
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
