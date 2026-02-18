#pragma once

#include <functional>

template <typename MfieldsState>
struct ExternalCurrentBase
{
  virtual ~ExternalCurrentBase() {}

  virtual void inject_current(MfieldsState& mflds) = 0;
};

template <typename MfieldsState>
struct ExternalCurrentFromLambda : ExternalCurrentBase<MfieldsState>
{
  ExternalCurrentFromLambda(std::function<void(MfieldsState&)> lambda)
    : lambda_{lambda}
  {}

  void inject_current(MfieldsState& mflds) override { return lambda_(mflds); }

private:
  std::function<void(MfieldsState&)> lambda_;
};
