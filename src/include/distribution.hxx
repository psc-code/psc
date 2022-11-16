#include <random>
#include <vector>
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

// ======================================================================
// InvertedCdf

template <typename Real>
struct InvertedCdf : public Distribution<Real>
{
  InvertedCdf(std::function<Real(Real)> cdf, int n_points = 100,
              Real eps = .0001)
  {
    Real xmin = -1, xmax = 1;

    while (cdf(xmin) < eps)
      xmin += .5;
    while (cdf(xmin) > eps)
      xmin -= 1 / 16.;
    while (cdf(xmax) > 1 - eps)
      xmax -= .5;
    while (cdf(xmax) < 1 - eps)
      xmax += 1 / 16.;

    x.push_back(xmin);
    this->cdf.push_back(0);
    Real dx = (xmax - xmin) / (n_points - 1);
    for (int i = 1; i < n_points - 1; i++) {
      x.push_back(xmin + i * dx);
      this->cdf.push_back(cdf(xmin + i * dx));
    }
    x.push_back(xmax);
    this->cdf.push_back(1);
  }

  Real get()
  {
    Real u = uniform.get();
    int guess_idx = cdf.size() * u;

    while (u < cdf[guess_idx]) {
      --guess_idx; // decrease cdf
    }
    while (u > cdf[guess_idx]) {
      ++guess_idx; // increase cdf
    }

    Real d_cdf = cdf[guess_idx] - cdf[guess_idx - 1];
    if (d_cdf == 0.)
      return x[guess_idx];

    // now: cdf[guess_idx-1] <= u < cdf[guess_idx]
    // so linearly interpolate between points

    Real w1 = (cdf[guess_idx] - u) / d_cdf;
    Real w2 = (u - cdf[guess_idx - 1]) / d_cdf;

    return w1 * x[guess_idx - 1] + w2 * x[guess_idx];
  }

private:
  std::vector<Real> cdf, x;
  Uniform<Real> uniform{0, 1};
};

} // namespace distribution