#include <random>

namespace distribution
{

// ======================================================================
// Distribution
// Base class for distributions

template <typename Real>
struct Distribution
{
  virtual Real get() = 0;
};

// ======================================================================
// Uniform

template <typename Real>
struct Uniform : public Distribution<Real>
{
  Uniform(Real min, Real max) : dist{min, max} {}
  Real get() override { return dist(generator); }

private:
  std::default_random_engine generator;
  std::uniform_real_distribution<Real> dist;
};

} // namespace distribution