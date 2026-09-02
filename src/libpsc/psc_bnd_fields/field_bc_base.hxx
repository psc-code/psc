#pragma once

template <typename MfieldsState>
struct FieldBcBase
{
  virtual ~FieldBcBase() {}

  virtual void apply_j_bcs(MfieldsState& mflds) = 0;
  virtual void apply_e_bcs(MfieldsState& mflds) = 0;
  virtual void apply_h_bcs(MfieldsState& mflds) = 0;
};
