
#include "gtest/gtest.h"

#include "testing.hxx"

// ======================================================================
// MomentTest

template <typename T>
struct MomentTest : ::testing::Test
{
  using Mfields = typename T::Mfields;
  using Mparticles = typename T::Mparticles;
  using dim_t = typename T::dim;
  using real_t = typename Mparticles::real_t;

  const double L = 160;

  const real_t fnqx = .05, fnqy = .05, fnqz = .05;
  const real_t dx = 10., dy = 10., dz = 10.;

  Int3 ibn = {2, 2, 2};

  void make_psc(const Grid_t::Kinds& kinds)
  {
    Int3 gdims = {16, 16, 16};
    if (dim_t::InvarX::value) {
      gdims[0] = 1;
      ibn[0] = 0;
    }
    if (dim_t::InvarY::value) {
      gdims[1] = 1;
      ibn[1] = 0;
    }
    if (dim_t::InvarZ::value) {
      gdims[2] = 1;
      ibn[2] = 0;
    }

    auto grid_domain = Grid_t::Domain{gdims, {L, L, L}};
    auto grid_bc =
      psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                    {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                    {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                    {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 200;
    auto coeff = Grid_t::Normalization{norm_params};

    grid_.reset(new Grid_t{grid_domain, grid_bc, kinds, coeff, 1., -1, ibn});
  }

  const Grid_t& grid()
  {
    assert(grid_);
    return *grid_;
  }

  std::unique_ptr<Grid_t> grid_;
};

template <typename MF, typename MP, typename D>
struct MomentTestConfig
{
  using Mfields = MF;
  using Mparticles = MP;
  using dim = D;
};

using MomentTestTypes = ::testing::Types<
#ifdef USE_CUDA
  MomentTestConfig<MfieldsCuda, MparticlesCuda<BS444>, dim_xyz>,
  MomentTestConfig<MfieldsCuda, MparticlesCuda<BS144>, dim_yz>,
#endif
  MomentTestConfig<MfieldsC, MparticlesDouble, dim_xyz>,
  MomentTestConfig<MfieldsC, MparticlesDouble, dim_yz>,
  MomentTestConfig<MfieldsSingle, MparticlesSingle, dim_xyz>,
  MomentTestConfig<MfieldsSingle, MparticlesSingle, dim_yz>>;

TYPED_TEST_SUITE(MomentTest, MomentTestTypes);

TYPED_TEST(MomentTest, Moment_n_1)
{
  using base_type = MomentTest<TypeParam>;
  using dim_t = typename base_type::dim_t;
  using Mfields = typename base_type::Mfields;
  using Mparticles = typename base_type::Mparticles;
  using real_t = typename base_type::real_t;

  using Moment = Moment_n_1st<Mfields, dim_t>;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {0., 0., 1.}, 1., 0});
  }
  Moment moment{mprts};
  auto gt = moment.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

template <typename Mfields>
struct MfieldsToHost
{
  using type = Mfields;
};

#ifdef USE_CUDA

template <>
struct MfieldsToHost<MfieldsCuda>
{
  using type = MfieldsSingle;
};

#endif

TYPED_TEST(MomentTest, Moments_1st)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using MfieldsHost = typename MfieldsToHost<Mfields>::type;
  using Moments = Moments_1st<Mparticles, MfieldsHost, dim_t>;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {0., 0., 1.}, 1., 0});
  }
  Moments moments{mprts};
  auto gt = moments.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_n_2) // FIXME, mostly copied
{
  using base_type = MomentTest<TypeParam>;
  using dim_t = typename base_type::dim_t;
  using Mfields = typename base_type::Mfields;
  using Mparticles = typename base_type::Mparticles;
  using real_t = typename base_type::real_t;

  using Moment = Moment_n_1st<Mfields, dim_t>;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{25., 5., 5.}, {0., 0., 1.}, 1., 0});
  }

  int i0 = 2;
  if (dim_t::InvarX::value)
    i0 = 0;

  Moment moment{mprts};
  auto gt = moment.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == i0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_v_1st)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Moments = Moment_v_1st<Mfields, dim_t>;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {.001, .002, .003}, 1., 0});
  }
  Moments moments{mprts};
  auto gt = moments.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, .005 * .001, eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_p_1st)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Moments = Moment_p_1st<Mfields, dim_t>;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {.001, .002, .003}, 1., 0});
  }
  Moments moments{mprts};
  auto gt = moments.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, .005 * .001, eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_rho_1st_nc_cc)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mfields::real_t;
  using Moment_t = Moment_rho_1st_nc<Mfields, dim_t>;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {0., 0., 0.}, 1., 0});
  }
  Moment_t moment{mprts};
  auto gt = moment.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if ((i == 0 && j == 0 && k == 0) || (i == 1 && j == 0 && k == 0) ||
            (i == 0 && j == 1 && k == 0) || (i == 1 && j == 1 && k == 0) ||
            (i == 0 && j == 0 && k == 1) || (i == 1 && j == 0 && k == 1) ||
            (i == 0 && j == 1 && k == 1) || (i == 1 && j == 1 && k == 1)) {
          EXPECT_NEAR(val, .005 / 8., eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if ((i == 0 && j == 0 && k == 0) || (i == 0 && j == 1 && k == 0) ||
            (i == 0 && j == 0 && k == 1) || (i == 0 && j == 1 && k == 1)) {
          EXPECT_NEAR(val, .005 / 4., eps)
            << "ijk " << i << " " << j << " " << k;

        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_rho_1st_nc_nc)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mfields::real_t;
  using Moment_t = Moment_rho_1st_nc<Mfields, dim_t>;

  const real_t eps = 1e-4;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{10., 10., 10.}, {0., 0., 0.}, 1., 0});
  }
  Moment_t moment{mprts};
  auto gt = moment.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if (i == 1 && j == 1 && k == 1) {
          EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if (j == 1 && k == 1) {
          EXPECT_NEAR(val, .005, eps) << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_n_2nd_nc)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Moment = Moment_n_2nd_nc<Mfields, dim_t>;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  const real_t w = .5f;
  const int nicell = 200; // FIXME, comes from testing.hxx
  const real_t cori = 1.f / nicell;

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {0., 0., 1.}, w, 0});
  }
  Moment moment{mprts};
  auto gt = moment.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if ((i == 0 || i == 1) && (j == 0 || j == 1) && (k == 0 || k == 1)) {
          EXPECT_NEAR(val, w * cori / 8., eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if ((j == 0 || j == 1) && (k == 0 || k == 1)) {
          EXPECT_NEAR(val, w * cori / 4., eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_rho_2nd_nc)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Moment = Moment_rho_2nd_nc<Mfields, dim_t>;
  using real_t = typename Mfields::real_t;

  const real_t eps = 1e-6;
  auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  const real_t w = .5f;
  const int nicell = 200; // FIXME, comes from testing.hxx
  const real_t cori = 1.f / nicell;

  // init particles
  Mparticles mprts{grid};
  {
    auto injector = mprts.injector();
    injector[0]({{5., 5., 5.}, {0., 0., 1.}, w, 0});
  }
  Moment moment{mprts};
  auto gt = moment.gt();
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if ((i == 0 || i == 1) && (j == 0 || j == 1) && (k == 0 || k == 1)) {
          EXPECT_NEAR(val, w * cori / 8., eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if ((j == 0 || j == 1) && (k == 0 || k == 1)) {
          EXPECT_NEAR(val, w * cori / 4., eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., eps) << "ijk " << i << " " << j << " " << k;
        }
      }
    });
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
