
#pragma once

#include "grid.hxx"

inline Grid_t make_grid()
{
  auto domain = Grid_t::Domain{{8, 4, 2},
			       {80.,  40., 20.}, {-40., -20., 0.},
			       {2, 2, 1}};
  std::vector<Int3> offs = {{0, 0, 0}, {4, 0, 0}};
  return Grid_t(domain, offs);
}

static Grid_t make_grid_1()
{
  auto domain = Grid_t::Domain{{8, 4, 2},
			       {10., 10., 10.}};
  std::vector<Int3> offs = { { 0, 0, 0 } };
  return Grid_t(domain, offs);
}

// ======================================================================
// MakeTestGrid
//
// Sample grid for testing, with 2 patches (incomplete coverage)

struct MakeTestGrid
{
  Grid_t operator()()
  {
    return make_grid();
  }
};

// ======================================================================
// MakeTestGrid1
//
// Sample grid for testing, with 1 patch only

struct MakeTestGrid1
{
  Grid_t operator()()
  {
    return make_grid_1();
  }
};

