#pragma once

#include <vector>
#include <cmath>
#include "stdlib.h"
#include <random>
#include <chrono>
#include <functional>
#include <memory>
#include "mpi.h"

namespace rng
{

namespace detail
{
using Generator = std::default_random_engine;
using GeneratorPtr = std::shared_ptr<Generator>;
GeneratorPtr process_generator;

int get_process_seed()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return 10000 * (rank + 1);
}

GeneratorPtr get_process_generator()
{
  if (!process_generator) {
    process_generator = std::make_shared<Generator>(get_process_seed());
  }
  return process_generator;
}

} // namespace detail

// ======================================================================
// Uniform

template <typename Real>
struct Uniform
{
  Uniform(Real min, Real max)
    : dist(min, max), gen(detail::get_process_generator())
  {}

  Uniform() : Uniform(0, 1) {}

  Real get() { return dist(*gen); }

private:
  detail::GeneratorPtr gen;
  std::uniform_real_distribution<Real> dist;
};

// ======================================================================
// Normal

template <typename Real>
struct Normal
{
  Normal(Real mean, Real stdev)
    : dist(mean, stdev), gen(detail::get_process_generator())
  {}
  Normal() : Normal(0, 1) {}

  Real get() { return dist(*gen); }
  // FIXME remove me, or make standalone func
  Real get(Real mean, Real stdev)
  {
    // should be possible to pass params to existing dist
    return std::normal_distribution<Real>(mean, stdev)(*gen);
  }

private:
  detail::GeneratorPtr gen;
  std::normal_distribution<Real> dist;
};

// ======================================================================
// InvertedCdf

template <typename Real>
void findPdfSupport(std::function<Real(Real)>& cdf, Real& xmin, Real& xmax,
                    Real tolerance, int max_n_iter)
{
  Real xmin_ub = 0, xmax_lb = 0;
  Real sf;
  int n_iter;

  // find where cfd > tolerance
  n_iter = 0;
  sf = 1;
  while (n_iter < max_n_iter && cdf(xmin_ub) < tolerance) {
    xmin_ub += (sf *= 1.5);
    n_iter++;
  }

  // find where cfd < 1-tolerance
  n_iter = 0;
  sf = 1;
  while (n_iter < max_n_iter && cdf(xmax_lb) > 1 - tolerance) {
    xmax_lb -= (sf *= 1.5);
    n_iter++;
  }

  Real xmin_lb = xmin_ub, xmax_ub = xmax_lb;

  // find where cdf < tolerance
  n_iter = 0;
  sf = 1;
  while (n_iter < max_n_iter && cdf(xmin_lb) > tolerance) {
    xmin_lb -= (sf *= 1.5);
    n_iter++;
  }

  // find where cdf > 1-tolerance
  n_iter = 0;
  sf = 1;
  while (n_iter < max_n_iter && cdf(xmax_ub) < 1 - tolerance) {
    xmax_ub += (sf *= 1.5);
    n_iter++;
  }

  // bisect search for where cdf = tolerance
  n_iter = 0;
  xmin = (xmin_lb + xmin_ub) / 2.;
  xmax = (xmax_lb + xmax_ub) / 2.;
  // logDiff measures the precision of the range relative to that of each bound
  Real logDiff =
    std::log2((xmax - xmin) / std::max(xmax_ub - xmax_lb, xmin_ub - xmin_lb));
  // if logDiff were recalculated, it would go up by ~1 with each iteration, so
  // logDiff(t) ~= n_iter(t) + logDiff(0)
  while (n_iter < max_n_iter && n_iter + logDiff < 10) {
    if (cdf(xmin) < tolerance) {
      xmin_lb = xmin;
    } else {
      xmin_ub = xmin;
    }
    if (cdf(xmax) > 1 - tolerance) {
      xmax_ub = xmax;
    } else {
      xmax_lb = xmax;
    }
    xmin = (xmin_lb + xmin_ub) / 2.;
    xmax = (xmax_lb + xmax_ub) / 2.;
    n_iter++;
  }
}

template <typename Real>
struct InvertedCdf
{
  InvertedCdf(std::function<Real(Real)> cdf, int n_points = 100,
              Real tolerance = .0001, int max_n_iter = 50)
  {
    Real xmin, xmax;
    findPdfSupport(cdf, xmin, xmax, tolerance, max_n_iter);

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

} // namespace rng