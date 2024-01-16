#pragma once

#include <cmath>
#include "params_parser.hxx"
#include "input_parser.hxx"

// ======================================================================
// getCalculatedBoxSize

// Calculates the radial distance where the spike in the exact distribution
// function ends, according to the equation (in paper units):
//    v_phi = max_v = 4*k*B*r^3 / (1 + 8*k*r^2)
// The exact solution can be decomposed into the difference of two Gaussians.
// The positive Gaussian has a mean of 0 and stdev of 1, independently of
// all parameters. "max_v" represents an upper limit for v_phi according
// to this term.
// The negative Gaussian has a mean given by the RHS of the equation. It drifts
// up, approaching a line with slope B/2. The negative Gaussian is the source
// of the spike.
double getSpikeSize(double B, double k)
{
  double max_v = 4.5;
  // solve cubic with linear coefficient = 0
  double a = 4. * k * B;
  double b = -8. * max_v * k;
  double d = -max_v;

  double p = -b / (3. * a);
  double t = -d / (2. * a);
  double q = p * p * p + t;
  double s = sqrt(t * (2. * q - t));

  double beta = .001; // FIXME don't hardcode this (see psc_bgk.cxx get_beta())
  double r = (std::cbrt(q + s) + std::cbrt(q - s) + p) * beta;
  return r;
}

// The hole size is determined empirically from the input data. This function
// finds where the electron density <= 1 + epsilon, where epsilon is small.
double getHoleSize(ParsedData& data)
{
  double epsilon = 1e-4;
  const int COL_RHO = 0;
  const int COL_NE = 1;
  for (int row = data.get_nrows() - 1; row >= 0; row--) {
    if (data[row][COL_NE] > 1 + epsilon)
      return data[row][COL_RHO];
  }
  throw "Unable to determine hole size.";
}

// Calculates a box size big enough to resolve the spike and the hole.
double getCalculatedBoxSize(double B, double k, ParsedData& data)
{
  double spike_size = getSpikeSize(B, k);
  double hole_size = getHoleSize(data);
  LOG_INFO("spike radius: %f\thole radius: %f\n", spike_size, hole_size);
  return 2 * std::max(spike_size, hole_size);
}

// ======================================================================
// PscBgkParams

struct PscBgkParams
{
  double box_size;     // physical length of region along y and z
                       // (overrides rel_box_size unless set to -1)
  double rel_box_size; // length of y and z dims in calculated units
  double Hx;           // strength of transverse magnetic field
  double q_i;          // ion charge
  double n_i;          // ion number density
  double m_i;          // ion mass
  double q_e;          // electron charge
  double m_e;          // electron mass
  int n_grid;          // number of grid cells
  int n_patches;       // number of patches
  int nicell;          // number of particles per gripdoint when density=1

  double k;  // a parameter for BGK solutions
  double h0; // a parameter for BGK solutions

  int fields_every;    // interval for pfd output
  int moments_every;   // interval for pfd_moments output
  int gauss_every;     // interval for gauss output/checking
  int particles_every; // interval for particle output

  // For Boltzmann ion cases
  bool do_ion; // whether or not to make ions follow B.D.
  double T_i;  // ion temperature

  // For modifying initial conditions
  double v_e_coef;     // multiplier for electron velocity
  double T_e_coef;     // multiplier for electron temperature
  bool reverse_v;      // whether or not to reverse electron velocity
  bool reverse_v_half; // whether or not to reverse electron velocity for y<0
  bool maxwellian;     // whether or not to use Maxwellian instead of exact sol

  // For 3D cases
  int n_grid_3;          // number of grid points in 3rd dimension
  double box_size_3;     // physical length of 3rd dimension
                         // (overrides rel_box_size_3 unless set to -1)
  double rel_box_size_3; // length of 3rd dimension in calculated units
  int n_patches_3;       // number of patches in 3rd dimension

  void loadParams(ParsedParams parsedParams, ParsedData& data)
  {
    box_size = parsedParams.getAndWarnOrDefault<double>("box_size", -1);
    rel_box_size = parsedParams.getOrDefault<double>("rel_box_size", 1);
    Hx = parsedParams.get<double>("H_x");
    q_i = parsedParams.get<double>("q_i");
    n_i = parsedParams.get<double>("n_i");
    m_i = parsedParams.get<double>("m_i");
    q_e = parsedParams.get<double>("q_e");
    m_e = parsedParams.get<double>("m_e");
    n_grid = parsedParams.get<int>("n_grid");
    n_patches = parsedParams.getAndWarnOrDefault<int>("n_patches", -1);
    if (n_patches <= 0)
      n_patches = n_grid / parsedParams.get<int>("n_cells_per_patch");
    nicell = parsedParams.getOrDefault<int>("nicell", 100);
    if (parsedParams.warnIfPresent("reverse_v", "Set v_e_coef to +-1 instead."))
      exit(0);
    reverse_v_half = parsedParams.get<bool>("reverse_v_half");

    fields_every = parsedParams.getOrDefault<int>("fields_every", 200);
    moments_every = parsedParams.getOrDefault<int>("moments_every", 200);
    gauss_every = parsedParams.getOrDefault<int>("gauss_every", 200);
    particles_every = parsedParams.getOrDefault<int>("particles_every", 0);

    do_ion = parsedParams.getOrDefault<bool>("ion", false);
    T_i = parsedParams.getOrDefault<double>("T_i", 0);

    v_e_coef = parsedParams.getOrDefault<double>("v_e_coef", 1);
    T_e_coef = parsedParams.getOrDefault<double>("T_e_coef", 1);

    maxwellian = parsedParams.getOrDefault<bool>("maxwellian", false);
    k = parsedParams.getOrDefault<double>("k", .1);
    h0 = parsedParams.getOrDefault<double>("h0", .9);

    n_grid_3 = parsedParams.getOrDefault<int>("n_grid_3", 1);
    box_size_3 = parsedParams.getAndWarnOrDefault<double>("box_size_3", -1);
    rel_box_size_3 = parsedParams.getOrDefault<double>("rel_box_size_3", 1);
    if (n_grid_3 < parsedParams.get<int>("n_cells_per_patch")) {
      n_patches_3 = 1;
    } else {
      n_patches_3 = parsedParams.getAndWarnOrDefault<int>("n_patches_3", -1);
      if (n_patches_3 <= 0)
        n_patches_3 = n_grid_3 / parsedParams.get<int>("n_cells_per_patch");
    }

    if (box_size <= 0)
      box_size = rel_box_size * getCalculatedBoxSize(Hx, k, data);
    if (box_size_3 <= 0)
      box_size_3 = rel_box_size_3 * getCalculatedBoxSize(Hx, k, data);
  }
};
