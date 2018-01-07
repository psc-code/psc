
#ifndef PSC_PUSH_PARTICLES_1VB_H
#define PSC_PUSH_PARTICLES_1VB_H

#include <psc_push_particles_private.h>

#include "fields.hxx"

template<typename MF, typename D>
struct push_p_config
{
  using mfields_t = MF;
  using dim_t = D;
};

template<typename C>
struct push_p_ops
{
  static void push_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base);
  static void stagger_mprts(struct psc_push_particles *push, struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base);
};

#endif
