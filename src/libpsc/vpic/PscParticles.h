
#pragma once

template<class ParticlesBase, class FA, class IA, class AA, class HA>
struct PscParticles : ParticlesBase
{
  typedef ParticlesBase Base;
  typedef PscParticles<Base, FA, IA, AA, HA> Particles;
  typedef FA FieldArray;
  typedef IA Interpolator;
  typedef AA Accumulator;
  typedef HA HydroArray;

  using Base::Base;
  using typename Base::iterator;
  using typename Base::const_iterator;
  using typename Base::Species;
  using typename Base::Particle;
  using typename Base::Grid;

};


