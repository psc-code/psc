
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

template<class Mparticles>
struct GetterCuda
{
  GetterCuda(Mparticles& mprts)
    : mprts_(mprts.template get_as<MparticlesSingle>()) {}

  typename MparticlesSingle::patch_t& operator[](int p) { return mprts_[p]; }

private:
  MparticlesSingle& mprts_;
};

template<class Mparticles,
	 typename = typename std::enable_if<!(std::is_same<Mparticles, MparticlesCuda<BS144>>::value ||
					      std::is_same<Mparticles, MparticlesCuda<BS444>>::value)>::type>
Getter<Mparticles> make_getter(Mparticles& mprts)
{
  return Getter<Mparticles>(mprts);
}

template<class Mparticles,
	 typename = typename std::enable_if<std::is_same<Mparticles, MparticlesCuda<BS144>>::value ||
					    std::is_same<Mparticles, MparticlesCuda<BS444>>::value>::type>
GetterCuda<Mparticles> make_getter(Mparticles& mprts)
{
  return GetterCuda<Mparticles>(mprts);
}

struct CurrentReference
{
  int m;
  Int3 pos;
  double val;
};

// ======================================================================
// class PushParticlesTest

template<typename T>
struct PushParticlesTest : ::testing::Test
{
  using dim = typename T::dim;
  using Mparticles = typename T::Mparticles;
  using Mfields = typename T::Mfields;
  using PushParticles = typename T::PushParticles;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mparticles::real_t;

  const real_t eps = 1e-5;
  const double L = 160;

  const real_t fnqx = .05, fnqy = .05, fnqz = .05;
  const real_t dx = 10., dy = 10., dz = 10.;
  
  Int3 ibn = { 2, 2, 2 };
  
  ~PushParticlesTest()
  {
    ppsc = NULL; // FIXME, should use psc_destroy(ppsc), or really, get rid of ppsc...
  }

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

  template<typename FUNC>
  void runSingleParticleTest(FUNC init_fields, particle_t prt0, particle_t prt1,
			     std::vector<CurrentReference> curr_ref = {})
  {
    auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
    make_psc(kinds);

    // init fields
    mflds = new Mfields{grid(), NR_FIELDS, ibn};
    SetupFields<Mfields>::set(*mflds, init_fields);

    // init particle
    mprts = new Mparticles{grid()};
    auto n_prts_by_patch = std::vector<uint>{1};
    SetupParticles<Mparticles>::setup_particles(*mprts, n_prts_by_patch, [&](int p, int n) -> typename Mparticles::particle_t {
	return prt0;
      });

    //mprts->dump("mprts.dump");
  
    // do one step
    PushParticles pushp_;
    pushp_.push_mprts(*mprts, *mflds);

    // check against reference
    for (auto& prt : make_getter(*mprts)[0]) {
      EXPECT_NEAR(prt.pxi, prt1.pxi, eps);
      EXPECT_NEAR(prt.pyi, prt1.pyi, eps);
      EXPECT_NEAR(prt.pzi, prt1.pzi, eps);
      EXPECT_NEAR(prt.qni_wni_, prt1.qni_wni_, eps);
      EXPECT_NEAR(prt.xi, prt1.xi, eps);
      EXPECT_NEAR(prt.yi, prt1.yi, eps);
      EXPECT_NEAR(prt.zi, prt1.zi, eps);
    }

    if (!curr_ref.empty()) {
      checkCurrent(curr_ref);
    }
  }

  void checkCurrent(std::vector<CurrentReference>& curr_ref)
  {
    auto mflds_ref = Mfields{grid(), NR_FIELDS, ibn};
    auto flds_ref = mflds_ref[0];
    for (auto& ref : curr_ref) {
      if (dim::InvarX::value) { ref.pos[0] = 0; }
      flds_ref(ref.m, ref.pos[0], ref.pos[1], ref.pos[2]) = ref.val;
    }

    auto flds = (*mflds)[0];
    this->grid().Foreach_3d(2, 2, [&](int i, int j, int k) {
	for (int m = JXI; m <= JZI; m++) {
	  auto val = flds(m, i,j,k);
	  auto val_ref = flds_ref(m, i,j,k);
	  EXPECT_NEAR(val, val_ref, eps) << "ijk " << i << " " << j << " " << k << " m " << m;
	}
      });
  }
  
  Vec3<double> push_x(const particle_t& prt0, particle_t& prt1)
  {
    Vec3<double> xi1 = { prt0.xi + vx(prt1),
			 prt0.yi + vy(prt1),
			 prt0.zi + vz(prt1) };
    
    if (!dim::InvarX::value) prt1.xi = xi1[0];
    if (!dim::InvarY::value) prt1.yi = xi1[1];
    if (!dim::InvarZ::value) prt1.zi = xi1[2];

    return xi1;
  }

  Mparticles* mprts;
  Mfields* mflds;
};

template<typename DIM, typename PUSHP, typename ORDER>
struct TestConfig
{
  using dim = DIM;
  using order = ORDER;
  using PushParticles = PUSHP;
  using Mparticles = typename PushParticles::Mparticles;
  using Mfields = typename PushParticles::Mfields;
  using Checks = Checks_<Mparticles, Mfields, ORDER>;
};

using TestConfig2ndDouble = TestConfig<dim_xyz,
				       PushParticles__<Config2ndDouble<dim_xyz>>,
				       checks_order_2nd>;
using TestConfig2ndDoubleYZ = TestConfig<dim_yz,
				       PushParticles__<Config2ndDouble<dim_yz>>,
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
using TestConfig1vbec3dCuda = TestConfig<dim_xyz,
					 PushParticlesCuda<CudaConfig1vbec3dGmem<dim_xyz, BS144>>,
					 checks_order_1st>;
using TestConfig1vbec3dCuda444 = TestConfig<dim_xyz,
					 PushParticlesCuda<CudaConfig1vbec3dGmem<dim_xyz, BS444>>,
					 checks_order_1st>;
using TestConfig1vbec3dCudaYZ = TestConfig<dim_yz,
					   PushParticlesCuda<CudaConfig1vbec3d<dim_yz, BS144>>,
					   checks_order_1st>;
#endif

using PushParticlesTestTypes = ::testing::Types<TestConfig2ndDoubleYZ,
						TestConfig1vbec3dSingleYZ,
#ifdef USE_CUDA
						TestConfig1vbec3dCudaYZ,
						TestConfig1vbec3dCuda,
						TestConfig1vbec3dCuda444,
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
// vx, vy, vz

template<typename particle_t>
typename particle_t::real_t vx(const particle_t& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
  return gamma * prt.pxi;
}

template<typename particle_t>
typename particle_t::real_t vy(const particle_t& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
  return gamma * prt.pyi;
}

template<typename particle_t>
typename particle_t::real_t vz(const particle_t& prt)
{
  auto gamma = 1./std::sqrt(1. + sqr(prt.pxi) + sqr(prt.pyi) + sqr(prt.pzi));
  return gamma * prt.pzi;
}

// ======================================================================
// SingleParticlePushp1 test
//
// zero EM test, check position push given existing pzi

TYPED_TEST(PushParticlesTest, SingleParticlePushp1)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.zi += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp2 test
//
// EZ = 1, check velocity push

TYPED_TEST(PushParticlesTest, SingleParticlePushp2)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return 2.;
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 3.;
  prt1.zi += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp3 test
//
// EZ = z, check interpolation in z direction
// 1vbec doesn't actually interpolate EZ in the z direction,
// but by choosing the midpoint, that works out
// (but this test isn't very comprehensive)

TYPED_TEST(PushParticlesTest, SingleParticlePushp3)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[2];
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 6.;
  prt1.zi += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp4 test
//
// EZ = y, check interpolation in y direction

TYPED_TEST(PushParticlesTest, SingleParticlePushp4)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[1];
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 5.; prt0.yi = 4.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 5.;
  prt1.zi = 5.980580;
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp5 test
//
// EZ = x, check interpolation in x direction

TYPED_TEST(PushParticlesTest, SingleParticlePushp5)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[0];
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 3.; prt0.yi = 5.; prt0.zi = 5.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 4.;
  if (Base::dim::InvarX::value) { prt1.pzi = 1.; }
  prt1.zi += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp6 test
//
// simply updating xi (in xyz, but not yz)

TYPED_TEST(PushParticlesTest, SingleParticlePushp6)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 1.; prt0.yi = 2.; prt0.zi = 3.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 1.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  if (!Base::dim::InvarX::value) prt1.xi += vx(prt1);
  prt1.yi += vy(prt1);
  prt1.zi += vz(prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp7 test
//
// test with EZ = z in different block (not just lower left)

TYPED_TEST(PushParticlesTest, SingleParticlePushp7)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    case EZ: return crd[2];
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 151.; prt0.yi = 152.; prt0.zi = 155.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 1.; prt0.pyi = 1.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  prt1.pzi = 156;
  this->push_x(prt0, prt1);
  
  this->runSingleParticleTest(init_fields, prt0, prt1);
}

// ======================================================================
// SingleParticlePushp8 test
//
// test current deposition in simple case

TYPED_TEST(PushParticlesTest, SingleParticlePushp8)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, this->fnqz / this->dz * (xi1[2] - prt0.zi) },
    };
  }
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp9 test
//
// test current deposition, z move only, but crossing bnd

TYPED_TEST(PushParticlesTest, SingleParticlePushp9)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 19.5;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 0.; prt0.pzi = 1.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JZI, {1, 1, 1}, this->fnqz / this->dz * (20. - prt0.zi) },
      { JZI, {1, 1, 2}, this->fnqz / this->dz * (xi1[2] - 20.) },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp10 test
//
// test current deposition, y move only, but crossing bnd

TYPED_TEST(PushParticlesTest, SingleParticlePushp10)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 19.5; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 0.; prt0.pyi = 1.; prt0.pzi = 0.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    auto fnqy = .05;
    auto dy = 10.;
    curr_ref = {
      { JYI, {1, 1, 1}, fnqy / dy * (20. - prt0.yi) },
      { JYI, {1, 2, 1}, fnqy / dy * (xi1[1] - 20.) },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
}

// ======================================================================
// SingleParticlePushp11 test
//
// test current deposition, x move only

TYPED_TEST(PushParticlesTest, SingleParticlePushp11)
{
  using Base = PushParticlesTest<TypeParam>;
  using particle_t = typename Base::particle_t;

  auto init_fields = [&](int m, double crd[3]) {
    switch (m) {
    default: return 0.;
    }
  };

  particle_t prt0, prt1;

  prt0.xi = 10.; prt0.yi = 10.; prt0.zi = 10.;
  prt0.qni_wni_ = 1.;
  prt0.pxi = 1.; prt0.pyi = 0.; prt0.pzi = 0.;
  prt0.kind_ = 0;

  prt1 = prt0;
  auto xi1 = this->push_x(prt0, prt1);

  std::vector<CurrentReference> curr_ref;
  if (std::is_same<typename TypeParam::order, checks_order_1st>::value) {
    curr_ref = {
      { JXI, {1, 1, 1}, this->fnqx / this->dx * (xi1[0] - prt0.xi) },
    };
  }
      
  this->runSingleParticleTest(init_fields, prt0, prt1, curr_ref);
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
  auto mflds = Mfields{grid, NR_FIELDS, this->ibn};
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
  auto mflds = Mfields{grid, NR_FIELDS, this->ibn};
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
