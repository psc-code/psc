#pragma once

#include "params_parser.hxx"

// ======================================================================
// PscBgkParams

struct PscBgkParams
{
  double box_size; // physical length of region along y and z
  double Hx;       // strength of transverse magnetic field
  double q_i;      // ion charge
  double n_i;      // ion number density
  double m_i;      // ion mass
  double q_e;      // electron charge
  double m_e;      // electron mass
  int n_grid;      // number of grid cells
  int n_patches;   // number of patches

  int fields_every;  // interval for pfd output
  int moments_every; // interval for pfd_moments output
  int gauss_every;   // interval for gauss output/checking

  // For Boltzmann ion cases
  bool do_ion; // whether or not to make ions follow B.D.
  double T_i;  // ion temperature

  // For modifying initial conditions
  double v_e_coef;     // multiplier for electron velocity
  double T_e_coef;     // multiplier for electron temperature
  bool reverse_v;      // whether or not to reverse electron velocity
  bool reverse_v_half; // whether or not to reverse electron velocity for y<0

  void loadParams(ParsedParams parsedParams)
  {
    box_size = parsedParams.get<double>("box_size");
    Hx = parsedParams.get<double>("H_x");
    q_i = parsedParams.get<double>("q_i");
    n_i = parsedParams.get<double>("n_i");
    m_i = parsedParams.get<double>("m_i");
    q_e = parsedParams.get<double>("q_e");
    m_e = parsedParams.get<double>("m_e");
    n_grid = parsedParams.get<int>("n_grid");
    n_patches = parsedParams.get<int>("n_patches");
    if (n_patches <= 0)
      n_patches = n_grid / parsedParams.get<int>("n_cells_per_patch");
    if (parsedParams.warnIfPresent("reverse_v", "Set v_e_coef to +-1 instead."))
      exit(0);
    reverse_v_half = parsedParams.get<bool>("reverse_v_half");

    fields_every = parsedParams.getOrDefault<int>("fields_every", 200);
    moments_every = parsedParams.getOrDefault<int>("moments_every", 200);
    gauss_every = parsedParams.getOrDefault<int>("gauss_every", 200);

    do_ion = parsedParams.getOrDefault<bool>("ion", false);
    T_i = parsedParams.getOrDefault<double>("T_i", 0);

    v_e_coef = parsedParams.getOrDefault<double>("v_e_coef", 1);
    T_e_coef = parsedParams.getOrDefault<double>("T_e_coef", 1);
  }
};