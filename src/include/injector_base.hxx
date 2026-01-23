#pragma once

template <typename MPARTICLES, typename MFIELDS_STATE>
struct InjectorBase
{
  using Mparticles = MPARTICLES;
  using MfieldsState = MFIELDS_STATE;

  virtual void inject(Mparticles& mprts, MfieldsState& mflds) = 0;
};

template <typename MPARTICLES, typename MFIELDS_STATE>
struct InjectFromLambda : InjectorBase<MPARTICLES, MFIELDS_STATE>
{
  using Mparticles = MPARTICLES;
  using MfieldsState = MFIELDS_STATE;

  InjectFromLambda(std::function<void(Mparticles&, MfieldsState&)> lambda)
    : lambda{lambda}
  {}

  void inject(Mparticles& mprts, MfieldsState& mflds) override
  {
    return lambda(mprts, mflds);
  }

private:
  std::function<void(Mparticles&, MfieldsState&)> lambda;
};
