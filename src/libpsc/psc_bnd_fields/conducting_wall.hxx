#pragma once

#include "psc.h"
#include "../axis.hxx"
#include "kg/Vec3.h"
#include "field_bc_base.hxx"
#include "../psc_bnd/psc_bnd_util.hxx"

namespace psc
{
namespace bnd
{
namespace field
{

template <typename Dim, typename MfieldsState>
struct ConductingWall : FieldBcBase<MfieldsState>
{
  using dim_t = Dim;

  ConductingWall(Axis d, LoHi lohi) : d{d}, lohi{lohi} {}

  void apply_j_bcs(MfieldsState& mflds) override
  {
    // todo
  }

  void apply_e_bcs(MfieldsState& mflds) override
  {
    // todo
  }

  void apply_h_bcs(MfieldsState& mflds) override
  {
    // todo
  }

  Axis d;
  LoHi lohi;
};

} // namespace field
} // namespace bnd
} // namespace psc
