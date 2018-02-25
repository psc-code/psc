
#ifndef PSC_PARTICLES_SINGLE_BY_KIND_H
#define PSC_PARTICLES_SINGLE_BY_KIND_H

#include "psc_particles_private.h"

#include "particles.hxx"

using particle_single_by_kind_real_t = float;

struct particle_single_by_kind_t
{
  using real_t = particle_single_by_kind_real_t;
};

struct psc_mparticles_single_by_kind : psc_mparticles_base
{
  using Base = psc_mparticles_base;
  using particle_t = particle_single_by_kind_t;

  using Base::Base;

  bk_mparticles *bkmprts;

  int get_n_prts() const override
  {
    return bk_mparticles_n_prts(bkmprts);
  }

  void get_size_all(uint *n_prts_by_patch) const override
  {
    bk_mparticles_size_all(bkmprts, n_prts_by_patch);
  }

  void reserve_all(const uint *n_prts_by_patch) override
  {
    bk_mparticles_reserve_all(bkmprts, n_prts_by_patch);
  }

  void resize_all(const uint *n_prts_by_patch) override
  {
    bk_mparticles_resize_all(bkmprts, n_prts_by_patch);
  }
};

#define psc_mparticles_single_by_kind(mprts)({				\
      mrc_to_subobj(mprts, struct psc_mparticles_single_by_kind);	\
})

// ======================================================================
// mparticles_single_by_kind_t

using mparticles_single_by_kind_t = PscMparticles<psc_mparticles_single_by_kind>;

#endif
