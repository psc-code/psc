#include <random>
#include <cmath>

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

// ======================================================================
// Normal

template <typename Real>
struct Normal : public Distribution<Real>
{
  Normal(Real mean, Real stdev) : mean{mean}, stdev{stdev} {}
  Normal() : Normal(0, 1) {}

  Real get() override { return get(mean, stdev); }
  Real get(Real mean, Real stdev)
  {
    if (has_next_value) {
      has_next_value = false;
      return next_value;
    }

    // Box-Muller transform to generate 2 independent values at once
    Real r = stdev * std::sqrt(-2 * std::log(1 - uniform.get()));
    Real theta = uniform.get() * 2 * M_PI;

    next_value = r * std::cos(theta) + mean;
    return r * std::sin(theta) + mean;
  }

private:
  Uniform<Real> uniform{0, 1};
  Real mean, stdev;
  Real next_value;
  bool has_next_value;
};

} // namespace distribution