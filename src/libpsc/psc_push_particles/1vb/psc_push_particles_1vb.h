
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include "push_particles.hxx"
#include "fields.hxx"

#include "../inc_defs.h"
#include "../push_config.hxx"

// ======================================================================
// PushParticles_

template<template<class> class PUSH_P_OPS>
class PushParticles_ : public PushParticlesBase
{
  using Self = PushParticles_<PUSH_P_OPS>;
  
public:
  PUSH_P_OPS<dim_yz> push_yz_;
};


#endif
