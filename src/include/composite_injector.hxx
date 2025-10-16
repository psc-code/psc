#pragma once

/// @brief An injector comprising two other injectors. Enables multiple
/// injectors in a single psc case.
/// @tparam INJECTOR_1 first type of injector
/// @tparam INJECTOR_2 second type of injector
template <typename INJECTOR_1, typename INJECTOR_2>
class CompositeInjector
{
public:
  using Mparticles = typename INJECTOR_1::Mparticles;
  using MfieldsState = typename INJECTOR_1::MfieldsState;

  CompositeInjector(INJECTOR_1 injector_1, INJECTOR_2 injector_2)
    : injector_1{injector_1}, injector_2{injector_2}
  {}

  /// @brief Calls both member injectors in order.
  /// @param mprts particles
  /// @param mflds fields
  void operator()(Mparticles& mprts, MfieldsState& mflds)
  {
    injector_1(mprts, mflds);
    injector_2(mprts, mflds);
  }

private:
  INJECTOR_1 injector_1;
  INJECTOR_2 injector_2;
};
