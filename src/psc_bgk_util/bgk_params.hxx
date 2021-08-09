#pragma once

#include "params_parser.hxx"

// ======================================================================
// PscBgkParams

struct PscBgkParams
{
  // physical length of region along y and z
  double box_size;
  // strength of transverse magnetic field
  double Hx;
  // ion charge
  double q_i;
  // ion number density
  double n_i;
  // ion mass
  double m_i;
  // electron charge
  double q_e;
  // electron mass
  double m_e;
  // number of grid cells
  int n_grid;
  // number of patches
  int n_patches;
  // whether or not to negate v (affects stability)
  bool reverse_v;
  // whether or not to negate v for y<0 (should destroy stability)
  bool reverse_v_half;

  // interval for pfd output
  int fields_every;
  // interval for pfd_moments output
  int moments_every;
  // interval for gauss output/checking
  int gauss_every;

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
    reverse_v = parsedParams.get<bool>("reverse_v");
    reverse_v_half = parsedParams.get<bool>("reverse_v_half");

    fields_every = parsedParams.getOrDefault<int>("fields_every", 200);
    moments_every = parsedParams.getOrDefault<int>("moments_every", 200);
    gauss_every = parsedParams.getOrDefault<int>("gauss_every", 200);
  }
};