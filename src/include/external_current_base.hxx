#pragma once

#include <functional>

template <typename MFIELDS_STATE>
struct ExternalCurrentBase
{
  using MfieldsState = MFIELDS_STATE;

  virtual ~ExternalCurrentBase() {}

  virtual void inject_current(MfieldsState& mflds) = 0;
};

template <typename MFIELDS_STATE>
struct ExternalCurrentFromLambda : ExternalCurrentBase<MFIELDS_STATE>
{
  using MfieldsState = MFIELDS_STATE;

  ExternalCurrentFromLambda(std::function<void(MfieldsState&)> lambda)
    : lambda_{lambda}
  {}

  void inject_current(MfieldsState& mflds) override { return lambda_(mflds); }

private:
  std::function<void(MfieldsState&)> lambda_;
};
