#pragma once

template <typename MfieldsState, typename Mparticles>
struct GaussCorrectorBase
{
  virtual void correct_gauss(MfieldsState& mflds, Mparticles& mprts) = 0;
};
