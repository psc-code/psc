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

} // namespace bnd
} // namespace psc
