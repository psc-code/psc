#pragma once

#include "grid.hxx"

namespace psc
{
namespace bnd
{

enum LoHi
{
  Lo,
  Hi,
};

/// whether patch `p` touches the domain's `lohi` boundary normal to `d`
inline bool at_boundary(const Grid_t& grid, int p, int d, LoHi lohi)
{
  return lohi == Lo ? grid.atBoundaryLo(p, d) : grid.atBoundaryHi(p, d);
}

/// whether particles of patch `p` can leave through the domain's `lohi`
/// boundary normal to `d`
inline bool at_open_particle_boundary(const Grid_t& grid, int p, int d,
                                      LoHi lohi)
{
  int bc = lohi == Lo ? grid.bc.prt_lo[d] : grid.bc.prt_hi[d];
  return bc == BND_PRT_OPEN && at_boundary(grid, p, d, lohi);
}

} // namespace bnd
} // namespace psc
