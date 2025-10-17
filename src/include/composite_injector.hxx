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

/// @brief Helper method to deduce the type of a composite injector, enabling
/// e.g. `auto composite_injector = make_composite(injector_1, injector_2);`
///
/// Wouldn't be necessary in C++17; see
/// https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rt-deduce.
/// @tparam INJECTOR_1 first type of injector
/// @tparam INJECTOR_2 second type of injector
/// @param injector_1 first injector
/// @param injector_2 second injector
/// @return
template <typename INJECTOR_1, typename INJECTOR_2>
CompositeInjector<INJECTOR_1, INJECTOR_2> make_composite(INJECTOR_1 injector_1,
                                                         INJECTOR_2 injector_2)
{
  return CompositeInjector<INJECTOR_1, INJECTOR_2>(injector_1, injector_2);
}
