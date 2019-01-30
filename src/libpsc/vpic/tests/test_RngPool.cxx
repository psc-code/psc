
#include "testing.h"

#include "VpicRng.h"
#include "PscRng.h"

#include "util/checkpt/checkpt.h"
#include "util/mp/mp.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>
#include <typeinfo>
#include "mrc_common.h"

// ----------------------------------------------------------------------
// test_Urng
//
// test that Rng behaves like a C++ Urng,
// though in the end we're not actually using this interface,
// at least for now

template<class G>
void test_Urng(G& g)
{
  std::cout << "-- test_Urng: " << typeid(g).name() << std::endl;
  g.seed(1);

  for (int i = 0; i < 3; i++) {
    std::cout << "gen " << i << " " << g() << std::endl;
  }
  
  for (int i = 0; i < 3; i++) {
    std::cout << "can   " << i << " " << std::generate_canonical<double, 2>(g) << std::endl;
  }

  std::uniform_real_distribution<double> uni(1, 2);
  for (int i = 0; i < 3; i++) {
    std::cout << "uni   " << i << " " << uni(g) << std::endl;
  }
  
}

template<class Rng>
void test_Urng()
{
  // std::mt19937 mt;
  // test_Urng(mt);

  Rng& rng = *(Rng::create());
  test_Urng(rng);
}

template<class RngPool, bool check_ref>
void test_RngPool()
{
  mprintf("-- test_RngPool: %s\n", typeid(RngPool).name());

  bool check = check_ref && (psc_world_rank == 0 ||
			     (psc_world_rank == 1 && psc_world_size == 2));
  
  RngPool pool;
  pool.seed(0, 0);
  auto rng = pool[0];

  // reference values for VPIC's Mersenne Twister and uniform real distribution
  double vals[2][3] = { {  0.93447864356331700186, -0.85443051718174234388, -0.41999871425083479259 },
			{ -0.44692117455683644245, -0.69362908232141173848, -0.66745219004201739033 }, };
  for (int i = 0; i < 3; i++) {
    double val = rng->uniform(-1., 1.);
    mprintf("uni %d %.20g\n", i, val);
    assert(val >= -1. && val <= 1.);
    if (check) {
      assert(val == vals[psc_world_rank][i]);
    }
  }

  // reference values for VPIC's Mersenne Twister and normal distribution
  double nvals[2][3] = { { 2.0612612026289509615, -1.7184880229355798953, 0.61942045083350327772 },
			 { 2.6605241085475155316, -1.22946085323322718  , 2.396498223337572675   }, };
  for (int i = 0; i < 3; i++) {
    double val = rng->normal(1., 2.);
    mprintf("nrm %d %.20g\n", i, val);
    if (check) {
      assert(val == nvals[psc_world_rank][i]);
    }
  }
}

int main(int argc, char **argv)
{
  testing_init(&argc, &argv);

  test_Urng<VpicRng>();
  test_Urng<PscRng>();

  // The VpicRng based pools should match the exact reference result sequence
  test_RngPool<VpicRngPool<VpicRng>, true>();
  test_RngPool<PscRngPool<VpicRng>, true>();
  // But PscRng, while still based on a Mersenne Twister, is using
  // C++-11 provided MT and distributions, so it won't match
  test_RngPool<PscRngPool<PscRng>, false>();

  MPI_Finalize();
  return 0;
}

