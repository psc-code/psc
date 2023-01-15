
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_single.h"
#include "psc_fields_single.h"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "../libpsc/psc_balance/psc_balance_impl.hxx"

#ifdef USE_CUDA
#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "psc_fields_cuda.h"
#include "cuda_base.hxx"
#endif

template <typename _Mparticles, typename _MfieldsState, typename _Mfields,
          typename _Balance = Balance_<_Mparticles, _MfieldsState, _Mfields>>
struct Config
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using Balance = _Balance;
};

using BalanceTestTypes = ::testing::Types<
  Config<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>,
  Config<MparticlesDouble, MfieldsStateDouble, MfieldsC>
#ifdef USE_CUDA
  ,
  Config<MparticlesCuda<BS144>, MfieldsStateCuda, MfieldsCuda,
         Balance_<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>>
#endif
  >;

TYPED_TEST_SUITE(BalanceTest, BalanceTestTypes);

// ======================================================================
// BalanceTest

template <typename T>
struct BalanceTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using Particle = typename Mparticles::Particle;

  BalanceTest()
  {
    auto domain =
      Grid_t::Domain{{1, 8, 16}, {10., 80., 160.}, {0., -40., -80.}, {1, 2, 2}};
    auto bc = psc::grid::BC{};
    auto kinds = Grid_t::Kinds{};
    auto norm = Grid_t::Normalization{};
    double dt = .1;
    grid_ = new Grid_t{domain, bc, kinds, norm, dt};

    grid_->kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  Mparticles mk_mprts()
  {
    Mparticles mprts(grid());
    mprts.define_species("test_species", 1., 1., 100, 10, 10, 0);
    return mprts;
  }

  template <typename _Mparticles>
  void inject_test_particles(_Mparticles& mprts,
                             const std::vector<uint>& n_prts_by_patch)
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      auto n_prts = n_prts_by_patch[p];
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

  const Grid_t& grid() { return *grid_; }

protected:
  Grid_t* grid_;
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(BalanceTest, Constructor)
{
  using Balance = typename TypeParam::Balance;
  auto balance = Balance{1., true};
}

#if 0
// ----------------------------------------------------------------------
// Initial1
//
// already balanced, nothing tbd

TYPED_TEST(BalanceTest, Initial1)
{
  using Balance = typename TypeParam::Balance;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto balance = Balance{0., true};

  auto n_prts_by_patch = std::vector<uint>{};
  if (size == 1) {
    n_prts_by_patch = {4, 4, 4, 4};
  } else if (size == 2) {
    n_prts_by_patch = {4, 4};
  } else {
    assert(0);
  }
  balance.initial(this->grid_, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// Initial2
//
// not balanced (on 2 procs)

TYPED_TEST(BalanceTest, Initial2)
{
  using Balance = typename TypeParam::Balance;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto balance = Balance{0., true};

  auto n_prts_by_patch = std::vector<uint>{};
  if (size == 1) {
    n_prts_by_patch = {4, 4, 4, 4};
  } else if (size == 2) {
    if (rank == 0) {
      n_prts_by_patch = {4, 4};
    } else {
      n_prts_by_patch = {4, 100};
    }
  } else {
    assert(0);
  }
  balance.initial(this->grid_, n_prts_by_patch);
}
#endif

// ----------------------------------------------------------------------
// Initial3
//
// not balanced (on 2 procs), actually rebalance some fields

TYPED_TEST(BalanceTest, Initial3)
{
  using MfieldsState = typename TypeParam::MfieldsState;
  using Mfields = typename TypeParam::Mfields;
  using Balance = typename TypeParam::Balance;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto balance = Balance{0., true};

  auto mflds_state = MfieldsState{this->grid()};
  auto mflds = Mfields{this->grid(), 3, {2, 2, 2}};

  auto n_prts_by_patch = std::vector<uint>{};
  if (size == 1) {
    n_prts_by_patch = {4, 4, 4, 4};
  } else if (size == 2) {
    if (rank == 0) {
      n_prts_by_patch = {4, 4};
    } else {
      n_prts_by_patch = {4, 100};
    }
  } else {
    assert(0);
  }
  balance.initial(this->grid_, n_prts_by_patch);
}

// ----------------------------------------------------------------------
// Every1
//
// not balanced (on 2 procs), so actually rebalances

TYPED_TEST(BalanceTest, Every1)
{
  using MfieldsState = typename TypeParam::MfieldsState;
  using Mfields = typename TypeParam::Mfields;
  using Balance = typename TypeParam::Balance;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto balance = Balance{0., true};

  auto mflds_state = MfieldsState{this->grid()};
  auto mflds = Mfields{this->grid(), 3, {2, 2, 2}};

  auto n_prts_by_patch = std::vector<uint>{};
  if (size == 1) {
    n_prts_by_patch = {4, 4, 4, 4};
  } else if (size == 2) {
    if (rank == 0) {
      n_prts_by_patch = {4, 4};
    } else {
      n_prts_by_patch = {4, 100};
    }
  } else {
    assert(0);
  }
  auto mprts = this->mk_mprts();
  this->inject_test_particles(mprts, n_prts_by_patch);

  balance(this->grid_, mprts);
}

TEST(Balance, best_mapping)
{
  const int n_procs = 3;
  const int n_patches = 10;
  std::vector<double> capability(n_procs, 1.);
  std::vector<double> loads(n_patches);
  std::iota(loads.begin(), loads.end(), 1.);

  std::cout << "capability ";
  std::copy(capability.begin(), capability.end(),
            std::ostream_iterator<double>(std::cout, " "));
  std::cout << "\n";

  std::cout << "loads ";
  std::copy(loads.begin(), loads.end(),
            std::ostream_iterator<double>(std::cout, " "));
  std::cout << "\n";

  auto n_patches_all = psc::balance::best_mapping({capability, loads});
  psc::balance::print_stats({capability, loads}, n_patches_all, true);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
#ifdef USE_CUDA
  cuda_base_init();
#endif

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
