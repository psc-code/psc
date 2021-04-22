
#include "gtest/gtest.h"

#include "testing.hxx"

#include "../libpsc/psc_push_fields/marder_impl.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/marder_cuda_impl.hxx"
#endif

template <typename T>
struct PushFieldsTest : PushParticlesTest<T>
{};

using PushFieldsTestTypes =
  ::testing::Types<TestConfig1vbec3dSingleYZ,
#ifdef USE_CUDA
                   TestConfig1vbec3dCudaYZ, TestConfig1vbec3dCuda444,
#endif
                   TestConfig1vbec3dSingle>;

TYPED_TEST_SUITE(PushFieldsTest, PushFieldsTestTypes);

// ======================================================================
// Pushf1

TYPED_TEST(PushFieldsTest, Pushf1)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using dim = typename TypeParam::dim;
  using PushFields = typename TypeParam::PushFields;

  const typename Mparticles::real_t eps = 1e-2;

  this->make_psc({});
  const auto& grid = this->grid();

  const double kz = 2. * M_PI / grid.domain.length[2];

  // init fields
  auto mflds = MfieldsState{grid};
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case EY: return sin(kz * crd[2]);
      default: return 0.;
    }
  });

  // run test
  PushFields pushf_;
  pushf_.push_H(mflds, 1., dim{});

  // check result
  auto&& h_gt = gt::host_mirror(mflds.gt());
  gt::copy(mflds.gt(), h_gt);
  auto bnd = mflds.ibn();
  auto&& i_gt =
    h_gt.view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]), _s(bnd[2], -bnd[2]));
  gt::launch<5, gt::space::host>(
    i_gt.shape(), [&](int i, int j, int k, int m, int p) {
      double z = grid.patches[p].z_cc(k);
      EXPECT_NEAR(i_gt(i, j, k, HX, p), kz * cos(kz * z), eps);
    });
}

TYPED_TEST(PushFieldsTest, Pushf2)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using dim = typename TypeParam::dim;
  using PushFields = typename TypeParam::PushFields;

  const typename Mparticles::real_t eps = 1e-2;

  this->make_psc({});
  const auto& grid = this->grid();

  const double ky = 2. * M_PI / grid.domain.length[1];

  // init fields
  auto mflds = MfieldsState{grid};
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case HX: return cos(ky * crd[1]);
      default: return 0.;
    }
  });

  // run test
  PushFields pushf_;
  pushf_.push_E(mflds, 1., dim{});

  // check result
  auto&& h_gt = gt::host_mirror(mflds.gt());
  gt::copy(mflds.gt(), h_gt);
  auto bnd = mflds.ibn();
  auto&& i_gt =
    h_gt.view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]), _s(bnd[2], -bnd[2]));
  gt::launch<5, gt::space::host>(
    i_gt.shape(), [&](int i, int j, int k, int m, int p) {
      double y = grid.patches[p].y_cc(j);
      EXPECT_NEAR(i_gt(i, j, k, EZ, p), ky * sin(ky * y), eps);
    });
}

// ======================================================================
// MarderCorrect

namespace detail
{
template <typename Mparticles, typename MfieldsState, typename Mfields,
          typename Dim>
struct marder_selector
{
  using type = Marder_<Mparticles, MfieldsState, Mfields>;
};

#ifdef USE_CUDA
template <typename Mparticles, typename Dim>
struct marder_selector<Mparticles, MfieldsStateCuda, MfieldsCuda, Dim>
{
  using type = MarderCuda<typename Mparticles::BS, Dim>;
};
#endif
} // namespace detail

template <typename Mparticles, typename MfieldsState, typename Mfields,
          typename Dim>
using marder_selector_t =
  typename detail::marder_selector<Mparticles, MfieldsState, Mfields,
                                   Dim>::type;

TYPED_TEST(PushFieldsTest, MarderCorrect)
{
  using Mparticles = typename TypeParam::Mparticles;
  using MfieldsState = typename TypeParam::MfieldsState;
  using Mfields = typename TypeParam::Mfields;
  using dim = typename TypeParam::dim;
  using PushFields = typename TypeParam::PushFields;

  const typename Mparticles::real_t eps = 1e-2;

  this->make_psc({});
  const auto& grid = this->grid();

  const double kz = 2. * M_PI / grid.domain.length[2];

  // run test
  double diffusion = .1;

  // init fields
  auto mflds = MfieldsState{grid};
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case EZ: return sin(kz * crd[2]);
      default: return 0.;
    }
  });
  auto mflds_ref = MfieldsState{grid};
  setupFields(mflds_ref, [&](int m, double crd[3]) {
    switch (m) {
      case EZ:
        return sin(kz * crd[2]) +
               .5 * grid.dt * diffusion * kz * cos(kz * crd[2]);
      default: return 0.;
    }
  });
  auto mphi = Mfields{grid, 1, grid.ibn};
  {
    auto&& h_mphi = hostMirror(mphi);
    double dz = grid.domain.dx[2];
    auto phi = make_Fields3d<dim_xyz>(h_mphi[0]);
    h_mphi.Foreach_3d(2, 2, [&](int i, int j, int k) {
      double z = k * dz;
      phi(0, i, j, k) = sin(kz * z);
    });
    copy(h_mphi, mphi);
  }

  psc::marder::correct(mflds, mphi, diffusion);

  // check result
  auto&& h_gt = gt::host_mirror(mflds.gt());
  auto&& h_gt_ref = gt::host_mirror(mflds_ref.gt());
  gt::copy(mflds.gt(), h_gt);
  gt::copy(mflds_ref.gt(), h_gt_ref);
  gt::launch<5, gt::space::host>(
    h_gt.shape(), [&](int i, int j, int k, int m, int p) {
      EXPECT_NEAR(h_gt(i, j, k, m, p), h_gt_ref(i, j, k, m, p), eps);
    });
}

template <typename T>
struct ItemTest : PushParticlesTest<T>
{};

using ItemTestTypes =
  ::testing::Types<TestConfig1vbec3dSingleYZ,
#ifdef USE_CUDA
                   TestConfig1vbec3dCudaYZ, TestConfig1vbec3dCuda444,
#endif
                   TestConfig1vbec3dSingle>;

TYPED_TEST_SUITE(ItemTest, ItemTestTypes);

// ======================================================================
// ItemDivE

TYPED_TEST(ItemTest, ItemDivE)
{
  using MfieldsState = typename TypeParam::MfieldsState;
  using Mfields = typename TypeParam::Mfields;
  using Item = Item_dive<MfieldsState>;

  const typename MfieldsState::real_t eps = 1e-2;

  this->make_psc({});
  const auto& grid = this->grid();

  double ky = 2. * M_PI / grid.domain.length[1];
  double kz = 2. * M_PI / grid.domain.length[2];

  // init fields
  auto mflds = MfieldsState{grid};
  setupFields(mflds, [&](int m, double crd[3]) {
    switch (m) {
      case EY: return cos(ky * crd[1]);
      case EZ: return sin(kz * crd[2]);
      default: return 0.;
    }
  });
  auto dx = grid.domain.dx;

  auto item_dive = Item(mflds);
  auto&& rho = gt::eval(item_dive.gt());

  auto rho_ref = gt::empty_like(rho);
  auto k_rho_ref = rho_ref.to_kernel();
  gt::launch<5, gt::space::host>(
    rho_ref.shape(), [=](int i, int j, int k, int m, int p) {
      double y = j * dx[1], z = k * dx[2];
      k_rho_ref(i, j, k, 0, p) = -ky * sin(ky * y) + kz * cos(kz * z);
    });

  // check result
  auto&& h_rho = host_mirror(rho);
  auto&& h_rho_ref = host_mirror(rho_ref);
  gt::copy(rho, h_rho);
  gt::copy(rho_ref, h_rho_ref);
  gt::launch<5, gt::space::host>(
    h_rho.shape(), [=](int i, int j, int k, int m, int p) {
      EXPECT_NEAR(h_rho(i, j, k, m, p), h_rho_ref(i, j, k, m, p), eps)
        << "i " << i << " j " << j << " k " << k << "\n";
    });
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
