
#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc_particles_single.h"
#include "psc_fields_single.h"
#include "psc_particles_double.h"
#include "psc_fields_c.h"
#include "../libpsc/psc_balance/psc_balance_impl.hxx"

template<typename _Mparticles, typename _MfieldsState, typename _Mfields>
struct Config
{
  using Mparticles = _Mparticles;
  using MfieldsState = _MfieldsState;
  using Mfields = _Mfields;
  using Balance = Balance_<Mparticles, MfieldsState, Mfields>;
};

using BalanceTestTypes = ::testing::Types<Config<MparticlesSingle, MfieldsStateSingle, MfieldsSingle>
					 ,Config<MparticlesDouble, MfieldsStateDouble, MfieldsC>
					 >;

TYPED_TEST_CASE(BalanceTest, BalanceTestTypes);

// ======================================================================
// BalanceTest

template<typename T>
struct BalanceTest : ::testing::Test
{
  using Mparticles = typename T::Mparticles;
  using Balance = typename T::Balance;
  using Particle = typename Mparticles::Particle;

  BalanceTest()
  {
    auto domain = Grid_t::Domain{{1, 8, 16},
				 {10., 80., 160.}, {0., -40., -80.},
				 {1, 2, 2}};
    auto bc = GridBc{};
    auto kinds = Grid_t::Kinds{};
    auto norm = Grid_t::Normalization{};
    double dt = .1;
    grid_ = new Grid_t{domain, bc, kinds, norm, dt};

    grid_->kinds.emplace_back(Grid_t::Kind(1., 1., "test_species"));
  }

  Mparticles mk_mprts()
  {
    Mparticles mprts(grid());
    mprts.define_species("test_species", 1., 1., 100, 10,
			 10, 0);
    return mprts;
  }

  template<typename _Mparticles>
  void inject_test_particles(_Mparticles& mprts, int n_prts)
  {
    auto inj = mprts.injector();
    for (int p = 0; p < mprts.n_patches(); ++p) {
      auto injector = inj[p];
      auto& patch = mprts.grid().patches[p];
      for (int n = 0; n < n_prts; n++) {
	double nn = double(n) / n_prts;
	auto L = patch.xe - patch.xb;
	particle_inject prt = {};
	prt.x[0] = patch.xb[0] + nn * L[0];
	prt.x[1] = patch.xb[1] + nn * L[1];
	prt.x[2] = patch.xb[2] + nn * L[2];
	prt.kind = 0;
	prt.w = 1.;
	injector(prt);
      }
    }
  }

  const Grid_t& grid() { return *grid_; }
  
protected:
  Grid_t *grid_;
};

// -----------------------------------------------------------------------
// Constructor

TYPED_TEST(BalanceTest, Constructor)
{
  using Balance = typename BalanceTest<TypeParam>::Balance;
  auto balance = Balance{1, 1., true};
}

// ----------------------------------------------------------------------
// Initial1
//
// already balanced, nothing tbd

TYPED_TEST(BalanceTest, Initial1)
{
  using Balance = typename BalanceTest<TypeParam>::Balance;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto balance = Balance{1, 0., true};

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
  using Balance = typename BalanceTest<TypeParam>::Balance;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  auto balance = Balance{1, 0., true};

  // auto mprts = this->mk_mprts();
  // this->inject_test_particles(mprts, n_prts);

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

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();

  MPI_Finalize();
  return rc;
}
