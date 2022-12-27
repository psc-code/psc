
#include "gtest/gtest.h"

#include "testing.hxx"

// ======================================================================
// Moment_rho_1st_nc_selector
//
// FIXME, should go away eventually

template <typename MF, typename D>
struct Moment_rho_1st_nc_selector
{
  using type = Moment_rho_1st_nc<typename MF::Storage, D>;
};

#ifdef USE_CUDA
template <typename D>
struct Moment_rho_1st_nc_selector<MfieldsCuda, D>
{
  using type = Moment_rho_1st_nc_cuda<D>;
};
#endif

// ======================================================================
// Moments_1st_selector
//
// FIXME, should go away eventually

template <typename MP, typename MF, typename D>
struct Moments_1st_selector
{
  using type = Moments_1st<MP, typename MF::Storage, D>;
};

#ifdef USE_CUDA
template <typename D>
struct Moments_1st_selector<MparticlesCuda<BS444>, MfieldsCuda, D>
{
  using type = Moments_1st_cuda<D>;
};
#endif

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
  const real_t eps = 1e-6;

  const real_t fnqx = .05, fnqy = .05, fnqz = .05;
  const real_t dx = 10., dy = 10., dz = 10.;
  const int nicell = 200;
  const real_t w = .4; // test particle weight

  Int3 ibn = {2, 2, 2};

  const Grid_t& make_grid()
  {
#ifdef USE_CUDA
    // if we're switching dim_yz <-> dim_xyz, cached maps become invalid
    BndCuda3::clear();
#endif

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
    norm_params.nicell = nicell;
    auto coeff = Grid_t::Normalization{norm_params};

    auto kinds = Grid_t::Kinds{Grid_t::Kind(1., 1., "test_species")};

    grid_.reset(new Grid_t{grid_domain, grid_bc, kinds, coeff, 1., -1, ibn});

    return *grid_;
  }

  Mparticles& make_mprts(const psc::particle::Inject& prt)
  {
    const auto& grid = make_grid();

    mprts_.reset(new Mparticles{grid});
    {
      auto injector = mprts_->injector();
      injector[0](prt);
    }
    return *mprts_;
  }

  const Grid_t& grid()
  {
    assert(grid_);
    return *grid_;
  }

  std::unique_ptr<Grid_t> grid_;
  std::unique_ptr<Mparticles> mprts_;
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
  using Moment = Moment_n_1st<typename Mfields::Storage, dim_t>;

  EXPECT_EQ(Moment::name(), "n_1st_cc");
  auto& mprts = this->make_mprts({{5., 5., 5.}, {0., 0., 1.}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  EXPECT_EQ(moment.comp_names(), std::vector<std::string>{"n_test_species"});
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, this->w / this->nicell, this->eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., this->eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moments_1st)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using MfieldsHost = hostMirror_t<Mfields>;
  using real_t = typename Mfields::real_t;
  using Moment =
    typename Moments_1st_selector<Mparticles, Mfields, dim_t>::type;

  EXPECT_EQ(Moment::name(), "all_1st_cc");
  auto& mprts = this->make_mprts({{5., 5., 5.}, {0., 0., 1.}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, this->w / this->nicell, this->eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., this->eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moments_1st_to_host)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using MfieldsHost = hostMirror_t<Mfields>;
  using real_t = typename Mfields::real_t;
  using Moment = Moments_1st<Mparticles, typename MfieldsHost::Storage, dim_t>;

  EXPECT_EQ(Moment::name(), "all_1st_cc");
  auto& mprts = this->make_mprts({{5., 5., 5.}, {0., 0., 1.}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, this->w / this->nicell, this->eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., this->eps) << "ijk " << i << " " << j << " " << k;
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
  using Moment = Moment_n_1st<typename Mfields::Storage, dim_t>;

  auto& mprts = this->make_mprts({{25., 5., 5.}, {0., 0., 1.}, this->w, 0});
  const auto& grid = this->grid();
  int i0 = 2;
  if (dim_t::InvarX::value)
    i0 = 0;

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == i0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, this->w / this->nicell, this->eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., this->eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_v_1st)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Moment = Moment_v_1st<typename Mfields::Storage, dim_t>;
  using real_t = typename Mfields::real_t;

  EXPECT_EQ(Moment::name(), "v_1st_cc");
  auto& mprts =
    this->make_mprts({{5., 5., 5.}, {.001, .002, .003}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, this->w / this->nicell * .001, this->eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., this->eps) << "ijk " << i << " " << j << " " << k;
      }
    });
  }
}

TYPED_TEST(MomentTest, Moment_p_1st)
{
  using Mparticles = typename TypeParam::Mparticles;
  using Mfields = typename TypeParam::Mfields;
  using dim_t = typename TypeParam::dim;
  using Moment = Moment_p_1st<typename Mfields::Storage, dim_t>;
  using real_t = typename Mfields::real_t;

  EXPECT_EQ(Moment::name(), "p_1st_cc");
  auto& mprts =
    this->make_mprts({{5., 5., 5.}, {.001, .002, .003}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (i == 0 && j == 0 && k == 0) {
        EXPECT_NEAR(val, this->w / this->nicell * .001, this->eps)
          << "ijk " << i << " " << j << " " << k;
      } else {
        EXPECT_NEAR(val, 0., this->eps) << "ijk " << i << " " << j << " " << k;
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
  using Moment = typename Moment_rho_1st_nc_selector<Mfields, dim_t>::type;

  EXPECT_EQ(Moment::name(), "rho_1st_nc");
  auto& mprts = this->make_mprts({{5., 5., 5.}, {0., 0., 0.}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if ((i == 0 && j == 0 && k == 0) || (i == 1 && j == 0 && k == 0) ||
            (i == 0 && j == 1 && k == 0) || (i == 1 && j == 1 && k == 0) ||
            (i == 0 && j == 0 && k == 1) || (i == 1 && j == 0 && k == 1) ||
            (i == 0 && j == 1 && k == 1) || (i == 1 && j == 1 && k == 1)) {
          EXPECT_NEAR(val, this->w / this->nicell / 8., this->eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., this->eps)
            << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if ((i == 0 && j == 0 && k == 0) || (i == 0 && j == 1 && k == 0) ||
            (i == 0 && j == 0 && k == 1) || (i == 0 && j == 1 && k == 1)) {
          EXPECT_NEAR(val, this->w / this->nicell / 4., this->eps)
            << "ijk " << i << " " << j << " " << k;

        } else {
          EXPECT_NEAR(val, 0., this->eps)
            << "ijk " << i << " " << j << " " << k;
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
  using Moment = typename Moment_rho_1st_nc_selector<Mfields, dim_t>::type;

  EXPECT_EQ(Moment::name(), "rho_1st_nc");
  auto& mprts = this->make_mprts({{10., 10., 10.}, {0., 0., 0.}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if (i == 1 && j == 1 && k == 1) {
          EXPECT_NEAR(val, this->w / this->nicell, this->eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., this->eps)
            << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if (j == 1 && k == 1) {
          EXPECT_NEAR(val, this->w / this->nicell, this->eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., this->eps)
            << "ijk " << i << " " << j << " " << k;
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
  using Moment = Moment_rho_2nd_nc<typename Mfields::Storage, dim_t>;
  using real_t = typename Mfields::real_t;

  EXPECT_EQ(Moment::name(), "rho_2nd_nc");
  auto& mprts = this->make_mprts({{5., 5., 5.}, {0., 0., 1.}, this->w, 0});
  const auto& grid = this->grid();

  Moment moment{grid};
  auto gt = psc::mflds::interior(grid, moment(mprts));
  for (int p = 0; p < grid.n_patches(); p++) {
    grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
      real_t val = gt(i, j, k, 0, p);
      if (std::is_same<dim_t, dim_xyz>::value) {
        if ((i == 0 || i == 1) && (j == 0 || j == 1) && (k == 0 || k == 1)) {
          EXPECT_NEAR(val, this->w / this->nicell / 8., this->eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., this->eps)
            << "ijk " << i << " " << j << " " << k;
        }
      } else if (std::is_same<dim_t, dim_yz>::value) {
        if ((j == 0 || j == 1) && (k == 0 || k == 1)) {
          EXPECT_NEAR(val, this->w / this->nicell / 4., this->eps)
            << "ijk " << i << " " << j << " " << k;
        } else {
          EXPECT_NEAR(val, 0., this->eps)
            << "ijk " << i << " " << j << " " << k;
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
