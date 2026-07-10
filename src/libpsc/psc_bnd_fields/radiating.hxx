#pragma once

#include "psc.h"
#include "../axis.hxx"
#include "kg/Vec3.h"
#include "kg/VecRange.hxx"
#include "field_bc_base.hxx"
#include "field_bc_util.hxx"
#include "../psc_bnd/psc_bnd_util.hxx"

namespace psc
{
namespace bnd
{
namespace field
{

template <typename Dim, typename MfieldsState, typename P>
struct Radiating : FieldBcBase<MfieldsState>
{
  using dim_t = Dim;
  using Pulse = P;

  Radiating(Pulse pulse, Axis d, LoHi lohi) : pulse{pulse}, d{d}, lohi{lohi} {}

  void apply_j_bcs(MfieldsState& mflds) override {}

  void apply_e_bcs(MfieldsState& mflds) override {}

  void apply_h_bcs(MfieldsState& mflds) override {}

  Axis d;
  LoHi lohi;

private:
  Pulse pulse;
};

} // namespace field
} // namespace bnd
} // namespace psc
