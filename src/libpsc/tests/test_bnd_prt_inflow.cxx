
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_double.h"
#include "psc_particles_single.h"
#include "particle_with_id.h"
#include "setup_particles.hxx"
#include "pushp.hxx"
#include "dim.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "../libpsc/cuda/mparticles_cuda.inl"
#endif
#include "particles_simple.inl"
#include <kg/io.h>

template <typename _Mparticles, typename _MakeGrid = MakeTestGrid1>
struct Config
{
  using Mparticles = _Mparticles;
  using MakeGrid = _MakeGrid;
};

using MparticlesInflowTestTypes =
  ::testing::Types<Config<MparticlesSingle>,
                   Config<MparticlesSingle, MakeTestGridYZ>
#ifdef USE_CUDA
                   ,
                   Config<MparticlesCuda<BS144>, MakeTestGridYZ1>,
                   Config<MparticlesCuda<BS144>, MakeTestGridYZ>
#endif
                   >;

TYPED_TEST_SUITE(MparticlesInflowTest, MparticlesInflowTestTypes);

// ======================================================================
// MparticlesInflowTest

template <typename T>
struct MparticlesInflowTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using Particle = typename Mparticles::Particle;
  using MakeGrid = typename T::MakeGrid;

  MparticlesInflowTest() : grid_{MakeGrid{}()}
  {
    grid_.kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  template <typename tag>
  Mparticles mk_mprts(tag dummy)
  {
    Mparticles mprts(grid_);
    mprts.define_species("test_species", 1., 1., 100, 10, 10, 0);
    return mprts;
  }

  Mparticles mk_mprts() { return mk_mprts(static_cast<Mparticles*>(nullptr)); }

  template <typename _Mparticles>
  void inject_test_particles(_Mparticles& mprts, int n_prts)
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
        double nn = double(n) / n_prts;
        auto L = patch.xe - patch.xb;
        psc::particle::Inject prt = {{patch.xb[0] + nn * L[0],
                                      patch.xb[1] + nn * L[1],
                                      patch.xb[2] + nn * L[2]},
                                     {},
                                     1.,
                                     0};
        injector(prt);
      }
    }
  }

  const Grid_t& grid() { return grid_; }

private:
  Grid_t grid_;
};

// ----------------------------------------------------------------------
// Inject

TYPED_TEST(MparticlesInflowTest, Inject)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject_test_particles(mprts, n_prts);
}

// -----------------------------------------------------------------------
// Inject2

TYPED_TEST(MparticlesInflowTest, Inject2)
{
  using Mparticles = typename TypeParam::Mparticles;
  const int n_prts = 1;

  auto mprts = this->mk_mprts();

  // FIXME, kinda the cheap way out
  if (mprts.n_patches() != 4) {
    return;
  }
  ASSERT_EQ(mprts.n_patches(), 4);

  int nn = 0;
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
        psc::particle::Inject prt{
          .5 * (patch.xb + patch.xe), {}, double(nn), 0};
        // use weight to store particle number for testing
        injector(prt);
        nn++;
      }
    }
  }

  EXPECT_EQ(mprts.size(), 4);

  auto n_prts_by_patch = mprts.sizeByPatch();
  EXPECT_EQ(n_prts_by_patch,
            std::vector<uint>({n_prts, n_prts, n_prts, n_prts}));

  nn = 0;
  auto accessor = mprts.accessor();
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto& patch = mprts.grid().patches[p];
    for (auto prt : accessor[p]) {
      // xm is patch-relative position
      auto xm = .5 * (patch.xe - patch.xb);
      EXPECT_EQ(prt.x()[0], xm[0]);
      EXPECT_EQ(prt.x()[1], xm[1]);
      EXPECT_EQ(prt.x()[2], xm[2]);
      EXPECT_EQ(prt.qni_wni(), nn);
      EXPECT_EQ(prt.w(), nn);
      EXPECT_EQ(prt.kind(), 0);

      auto x = .5 * (patch.xb + patch.xe);
      EXPECT_EQ(prt.position()[0], x[0]);
      EXPECT_EQ(prt.position()[1], x[1]);
      EXPECT_EQ(prt.position()[2], x[2]);
      nn++;
    }
  }
}

// ----------------------------------------------------------------------
// Injector

TYPED_TEST(MparticlesInflowTest, Injector)
{
  const int n_prts = 4;

  auto mprts = this->mk_mprts();

  this->inject_test_particles(mprts, n_prts);

  auto accessor = mprts.accessor();
  for (int p = 0; p < mprts.n_patches(); ++p) {
    auto prts = accessor[p];
    auto& patch = mprts.grid().patches[p];
    EXPECT_EQ(prts.size(), n_prts);
    int n = 0;
    for (auto prt : prts) {
      double nn = double(n) / n_prts;
      auto L = patch.xe - patch.xb;
      auto x = prt.position();
      EXPECT_EQ(x[0], patch.xb[0] + nn * L[0]);
      EXPECT_EQ(x[1], patch.xb[1] + nn * L[1]);
      EXPECT_EQ(x[2], patch.xb[2] + nn * L[2]);
      EXPECT_EQ(prt.w(), 1.);
      EXPECT_EQ(prt.kind(), 0);
      ++n;
    }
  }
}

TEST(TestSetupParticlesInflow, Simple)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, .1};
  Mparticles mprts{grid};

  SetupParticles<Mparticles> setup_particles(grid);

  Double3 pos = {5., 5., 5.};
  Double3 u = Double3{0.0, 0.0, 0.0};
  int tag = 0;
  psc_particle_np np = {0, 1.0, [u]() { return u; }, tag};
  double wni = 1.0;
  auto prt = setup_particles.setupParticle(np, pos, wni);

  EXPECT_EQ(prt.x, pos);
  EXPECT_EQ(prt.u, u);
  EXPECT_EQ(prt.w, wni);
  EXPECT_EQ(prt.kind, 0);
  EXPECT_EQ(prt.tag, tag);
}

TEST(TestSetupParticlesInflow, Maxwellian)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, .1};
  Mparticles mprts{grid};

  SetupParticles<Mparticles> setup_particles(grid);

  Double3 pos = {5., 5., 5.};
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  int tag = 0;

  psc_particle_npt npt = {0, 1.0, u, T, tag};
  auto p = setup_particles.createMaxwellian(npt);
  psc_particle_np np = {0, 1.0, p, tag};
  double wni = 1.0;
  auto prt = setup_particles.setupParticle(np, pos, wni);

  EXPECT_EQ(prt.x, pos);
  EXPECT_EQ(prt.u, u);
  EXPECT_EQ(prt.w, wni);
  EXPECT_EQ(prt.kind, 0);
  EXPECT_EQ(prt.tag, tag);
}

template <typename MP>
class Inflow
{
public:
  using Mparticles = MP;
  using real_t = typename Mparticles::real_t;

  // offset_in_cell: () -> double
  Inflow(const Grid_t& grid, psc_particle_npt npt,
         double (*offset_in_cell_dist)(), int inj_dim_idx)
    : grid_(grid),
      setup_particles_(grid),
      advance_(grid.dt),
      // FIXME np.p is a std::function and is called many times; better to use a
      // lambda
      np_{npt.kind, npt.n, setup_particles_.createMaxwellian(npt)},
      offset_in_cell_dist_(offset_in_cell_dist),
      inj_dim_idx_(inj_dim_idx)
  {}

  auto get_advanced_prt(Double3 pos, real_t wni)
  {
    auto prt = setup_particles_.setupParticle(np_, pos, wni);

    auto v = advance_.calc_v(prt.u);
    advance_.push_x(prt.x, v, 1.0);

    return prt;
  }

  template <typename Injector>
  void inject_into_boundary_cell(Injector& injector,
                                 Int3 boundary_cell_global_idx)
  {
    assert(boundary_cell_global_idx[inj_dim_idx_] == 0);
    boundary_cell_global_idx[inj_dim_idx_] = -1;

    int n_in_cell = setup_particles_.get_n_in_cell(np_);
    double wni = setup_particles_.getWeight(np_, n_in_cell);

    for (int cnt = 0; cnt < n_in_cell; cnt++) {
      Double3 offset = {offset_in_cell_dist_(), offset_in_cell_dist_(),
                        offset_in_cell_dist_()};
      auto pos =
        (Double3(boundary_cell_global_idx) + offset) * grid_.domain.dx +
        grid_.domain.corner;
      auto prt = get_advanced_prt(pos, wni);

      if (prt.x[inj_dim_idx_] < grid_.domain.corner[inj_dim_idx_]) {
        continue;
      }

      injector(prt);
    }
  }

  template <typename Injector>
  void inject_into_boundary_patch(Injector& injector,
                                  const Grid_t::Patch& boundary_patch)
  {
    Int3 ilo = boundary_patch.off;
    Int3 ihi = ilo + grid_.ldims;

    assert(ilo[inj_dim_idx_] == 0);

    int dim1 = std::min((inj_dim_idx_ + 1) % 3, (inj_dim_idx_ + 2) % 3);
    int dim2 = std::max((inj_dim_idx_ + 1) % 3, (inj_dim_idx_ + 2) % 3);

    for (Int3 cell_idx = ilo; cell_idx[dim1] < ihi[dim1]; cell_idx[dim1]++) {
      for (cell_idx[dim2] = ilo[dim2]; cell_idx[dim2] < ihi[dim2];
           cell_idx[dim2]++) {
        inject_into_boundary_cell(injector, cell_idx);
      }
    }
  }

  template <typename Injectors>
  void inject_into_boundary(Injectors& injectors_by_patch)
  {
    for (int patch_idx = 0; patch_idx < grid_.n_patches(); patch_idx++) {
      const auto& patch = grid_.patches[patch_idx];

      if (patch.off[inj_dim_idx_] != 0) {
        continue;
      }

      auto& injector = injectors_by_patch[patch_idx];
      inject_into_boundary_patch(injector, patch);
    }
  }

  const Grid_t& grid_;
  AdvanceParticle<real_t, dim_y> advance_;
  SetupParticles<Mparticles> setup_particles_;
  psc_particle_np np_;
  double (*offset_in_cell_dist_)();
  int inj_dim_idx_;
};

class TestInjector
{
public:
  void operator()(psc::particle::Inject prt) { prts.push_back(prt); }

  std::vector<psc::particle::Inject> prts;
};

double half() { return 0.5; }

TEST(TestSetupParticlesInflow, Advance)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 pos = {5., -2., 5.};
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles> inflow(grid, npt, *half, 1);

  auto prt = inflow.get_advanced_prt(pos, 1.0);

  EXPECT_NEAR(prt.x[0], 5., 1e-5);
  EXPECT_NEAR(prt.x[1], 2.47214, 1e-5);
  EXPECT_NEAR(prt.x[2], 5., 1e-5);
}

TEST(TestSetupParticlesInflow, InjectIntoCell)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.8, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles> inflow(grid, npt, *half, 1);

  TestInjector injector;
  inflow.inject_into_boundary_cell(injector, {0, 0, 0});

  EXPECT_EQ(injector.prts.size(), prm.nicell);

  for (auto prt : injector.prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);
    EXPECT_NEAR(prt.x[2], 5., 1e-5);
  }
}

TEST(TestSetupParticlesInflow, InjectIntoCellFilter)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  // too slow to enter domain
  Double3 u = {0.0, 0.5, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles> inflow(grid, npt, *half, 1);

  TestInjector injector;
  inflow.inject_into_boundary_cell(injector, {0, 0, 0});

  EXPECT_EQ(injector.prts.size(), 0);
}

TEST(TestSetupParticlesInflow, InjectIntoPatch)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 1, 1}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.8, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles> inflow(grid, npt, *half, 1);

  TestInjector injector;
  inflow.inject_into_boundary_patch(injector, grid.patches[0]);

  EXPECT_EQ(injector.prts.size(),
            domain.ldims[0] * domain.ldims[2] * prm.nicell);

  for (int i = 0; i < injector.prts.size(); i++) {
    auto prt = injector.prts[i];

    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);

    if (i < prm.nicell) {
      EXPECT_NEAR(prt.x[2], 5., 1e-5);
    } else {
      EXPECT_NEAR(prt.x[2], 15., 1e-5);
    }
  }
}
TEST(TestSetupParticlesInflow, InjectIntoBoundary)
{
  using Mparticles = MparticlesDouble;

  auto domain = Grid_t::Domain{{1, 2, 2}, {10., 20., 20.}, {}, {1, 2, 2}};
  double dt = 10.;
  auto kinds = Grid_t::Kinds{{1., 100., "i"}};
  auto prm = Grid_t::NormalizationParams::dimensionless();
  prm.nicell = 2;
  Grid_t grid{domain, {}, kinds, {prm}, dt};
  Mparticles mprts{grid};

  Double3 u = {0.0, 0.8, 0.0};
  Double3 T = {0.0, 0.0, 0.0};
  psc_particle_npt npt = {0, 1.0, u, T};

  Inflow<Mparticles> inflow(grid, npt, *half, 1);

  std::vector<TestInjector> injectors = {TestInjector(), TestInjector(),
                                         TestInjector(), TestInjector()};
  inflow.inject_into_boundary(injectors);

  EXPECT_EQ(injectors[0].prts.size(),
            domain.ldims[0] * domain.ldims[2] * prm.nicell);
  for (auto prt : injectors[0].prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);
    EXPECT_NEAR(prt.x[2], 5., 1e-5);
  }

  EXPECT_EQ(injectors[1].prts.size(), 0);

  EXPECT_EQ(injectors[2].prts.size(),
            domain.ldims[0] * domain.ldims[2] * prm.nicell);
  for (auto prt : injectors[2].prts) {
    EXPECT_NEAR(prt.x[0], 5., 1e-5);
    EXPECT_NEAR(prt.x[1], 1.246950, 1e-5);
    EXPECT_NEAR(prt.x[2], 15., 1e-5);
  }

  EXPECT_EQ(injectors[3].prts.size(), 0);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
