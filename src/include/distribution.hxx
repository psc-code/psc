#include <vector>
#include <cmath>
#include "stdlib.h"
#include <random>
#include <chrono>

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
  Uniform(Real min, Real max)
    : dist(min, max),
      gen(std::chrono::system_clock::now().time_since_epoch().count())
  {}

  Real get() override { return dist(gen); }

private:
  std::default_random_engine gen;
  std::uniform_real_distribution<Real> dist;
};

// ======================================================================
// Normal

template <typename Real>
struct Normal : public Distribution<Real>
{
  Normal(Real mean, Real stdev)
    : dist(mean, stdev),
      gen(std::chrono::system_clock::now().time_since_epoch().count())
  {}
  Normal() : Normal(0, 1) {}

  Real get() override { return dist(gen); }
  // FIXME remove me, or make standalone func
  Real get(Real mean, Real stdev)
  {
    // should be possible to pass params to existing dist
    return std::normal_distribution<Real>(mean, stdev)(gen);
  }

private:
  std::default_random_engine gen;
  std::normal_distribution<Real> dist;
};

// ======================================================================
// InvertedCdf

template <typename Real>
struct InvertedCdf : public Distribution<Real>
{
  InvertedCdf(std::function<Real(Real)> cdf, int n_points = 100,
              Real eps = .0001, int max_n_iter = 50)
  {
    Real xmin = -1, xmax = 1;

    for (int n_iter = 0; cdf(xmin) < eps && n_iter < max_n_iter; ++n_iter)
      xmin += .5;
    for (int n_iter = 0; cdf(xmin) > eps && n_iter < max_n_iter; ++n_iter)
      xmin -= 1 / 16.;
    for (int n_iter = 0; cdf(xmax) > 1 - eps && n_iter < max_n_iter; ++n_iter)
      xmax -= .5;
    for (int n_iter = 0; cdf(xmax) < 1 - eps && n_iter < max_n_iter; ++n_iter)
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