#pragma once

#include <functional>

template <typename Mparticles, typename MfieldsState>
struct DiagnosticBase
{
  virtual ~DiagnosticBase() {}

  virtual void perform_diagnostic(Mparticles& mprts, MfieldsState& mflds) = 0;
};

template <typename Mparticles>
struct ParticleDiagnosticBase
{
  virtual ~ParticleDiagnosticBase() {}

  virtual void perform_diagnostic(Mparticles& mprts) = 0;
};

template <typename Mparticles, typename MfieldsState>
struct DiagnosticFromLambda : public DiagnosticBase<Mparticles, MfieldsState>
{
  DiagnosticFromLambda(std::function<void(Mparticles&, MfieldsState&)> lambda)
    : lambda_{lambda}
  {}

  void perform_diagnostic(Mparticles& mprts, MfieldsState& mflds) override
  {
    return lambda_(mprts, mflds);
  }

private:
  std::function<void(Mparticles&, MfieldsState&)> lambda_;
};
