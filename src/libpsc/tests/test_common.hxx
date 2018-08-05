
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

