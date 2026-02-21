#pragma once

#include <functional>

template <typename Mparticles, typename MfieldsState>
struct InjectorBase
{
  virtual ~InjectorBase() {}

  virtual void inject(Mparticles& mprts, MfieldsState& mflds) = 0;
};

template <typename Mparticles, typename MfieldsState>
struct InjectFromLambda : InjectorBase<Mparticles, MfieldsState>
{
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
