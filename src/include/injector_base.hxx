#pragma once

#include <functional>

template <typename MPARTICLES, typename MFIELDS_STATE>
struct InjectorBase
{
  using Mparticles = MPARTICLES;
  using MfieldsState = MFIELDS_STATE;

  virtual ~InjectorBase() {}

  virtual void inject(Mparticles& mprts, MfieldsState& mflds) = 0;
};

template <typename MPARTICLES, typename MFIELDS_STATE>
struct InjectFromLambda : InjectorBase<MPARTICLES, MFIELDS_STATE>
{
  using Mparticles = MPARTICLES;
  using MfieldsState = MFIELDS_STATE;

  InjectFromLambda(std::function<void(Mparticles&, MfieldsState&)> lambda)
    : lambda_{lambda}
  {}

  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    return lambda_(mprts, mflds);
  }

private:
  std::function<void(Mparticles&, MfieldsState&)> lambda_;
};
