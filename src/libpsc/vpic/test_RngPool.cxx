
#include "Rng.h"

#include "util/checkpt/checkpt.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>


// ----------------------------------------------------------------------
// test_Rng_1
//
// test that Rng behaves like a C++ urng,
// though in the end we're not actually using this interface,
// at least for now

template<class G>
void test_Rng_1(G& g)
{
  std::cout << "-- test_Rng_1: " << typeid(g).name() << std::endl;
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
void test_Rng()
{
  std::mt19937 mt;
  test_Rng_1(mt);

  Rng& rng = *(Rng::create());
  test_Rng_1(rng);
}

template<class RngPool>
void test_RngPool()
{
  std::cout << "-- test_RngPool: " << typeid(RngPool).name() << std::endl;
  
  RngPool pool;
  auto rng = pool[0];

  // reference values for VPIC's Mersenne Twister and uniform real distribution
  double vals[3] = { 0.93447864356331700186, -0.85443051718174234388, -0.41999871425083479259 };
  for (int i = 0; i < 3; i++) {
    double val = rng->uniform(-1., 1.);
    std::cout << "uni " << i << " " << std::setprecision(20) << val << std::endl;
    assert(val >= -1. && val <= 1.);
    assert(val == vals[i]);
  }

  // reference values for VPIC's Mersenne Twister and normal distribution
  double nvals[3] = { 2.0612612026289509615, -1.7184880229355798953, 0.61942045083350327772 };
  for (int i = 0; i < 3; i++) {
    double val = rng->normal(1., 2.);
    std::cout << "nrm " << i << " " << std::setprecision(20) << val << std::endl;
    assert(val == nvals[i]);
  }
}

int main(int argc, char **argv)
{
  boot_checkpt(&argc, &argv);
  
  test_Rng<VpicRng>();
  
  test_RngPool<VpicRngPool<VpicRng>>();
}

