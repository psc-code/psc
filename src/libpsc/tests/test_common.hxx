
#pragma once

#include "grid.hxx"

// ======================================================================
// MakeTestGrid1
//
// Sample grid for testing, with 1 patch only

struct MakeTestGrid1
{
  Grid_t operator()()
  {
    auto domain = Grid_t::Domain{{8, 4, 2},
				 {10., 10., 10.}};
    std::vector<Int3> offs = { { 0, 0, 0 } };
    return Grid_t{domain, offs};
  }
};

// ======================================================================
// MakeTestGridYZ
//
// Sample grid for testing, with 2 patches (incomplete coverage)
// patch size 4 x 8 (2 4x4 blocks)
// p1 : [0:10] x [-40:0] x [-80:0]
// p2 : [0:10] x [-40:0] x [0:80]

struct MakeTestGridYZ
{
  Grid_t operator()()
  {
    auto domain = Grid_t::Domain{{1, 8, 16},
				 {10., 80., 160.}, {0., -40., -80.},
				 {1, 2, 2}};
    std::vector<Int3> offs = {{0, 0, 0}, {0, 0, 4}};
    return Grid_t{domain, offs};
  }
};

// ======================================================================
// MakeTestGridYZ1
//
// Sample grid for testing, with 1 patch
// patch size 8 x 16 (8 4x4 blocks)

struct MakeTestGridYZ1
{
  Grid_t operator()()
  {
    auto domain = Grid_t::Domain{{1, 8, 16},
				 {10., 80., 160.}, {0., -40., -80.},
				 {1, 1, 1}};
    std::vector<Int3> offs = {{0, 0, 0}};
    return Grid_t{domain, offs};
  }
};

