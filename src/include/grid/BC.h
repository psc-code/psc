
#ifndef GRID_BC_H
#define GRID_BC_H

// FIXME, enums should be namespaced, enum class

/// Possible boundary conditions for fields
enum
{
  BND_FLD_OPEN,
  BND_FLD_PERIODIC,
  BND_FLD_CONDUCTING_WALL,
  BND_FLD_ABSORBING,
};

/// Possible boundary conditions for particles
enum
{
  BND_PRT_REFLECTING,
  BND_PRT_PERIODIC,
  BND_PRT_ABSORBING,
  BND_PRT_OPEN,
};

namespace psc
{
namespace grid
{

/// Describes the spatial domain to operate on.
///
/// This struct describes the spatial dimension of the simulation-box
///@note Here, you can also set the dimensionality by eliminating a dimension.
/// Example: To simulate in xy only, set \verbatim psc_domain.gdims[2]=1
///\endverbatim Also, set the boundary conditions for the eliminated dimensions
/// to BND_FLD_PERIODIC or you'll get invalid \a dt and \a dx

struct BC
{
  BC()
    : fld_lo{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
      fld_hi{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
      prt_lo{BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
      prt_hi{BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}
  {}

  BC(Int3 fld_lo, Int3 fld_hi, Int3 prt_lo, Int3 prt_hi)
    : fld_lo(fld_lo), fld_hi(fld_hi), prt_lo(prt_lo), prt_hi(prt_hi)
  {}

  Int3
    fld_lo; ///< Boundary conditions of the fields. Can be any value of BND_FLD.
  Int3
    fld_hi; ///< Boundary conditions of the fields. Can be any value of BND_FLD.
  Int3 prt_lo; ///< Boundary conditions of the particles. Can be any value of
               ///< BND_PART.
  Int3 prt_hi; ///< Boundary conditions of the particles. Can be any value of
               ///< BND_PART.
};

} // namespace grid
} // namespace psc

#endif
