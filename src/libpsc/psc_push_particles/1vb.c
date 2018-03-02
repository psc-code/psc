
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "1vb/psc_push_particles_1vb.h"

#include "fields.hxx"

using real_t = mparticles_t::real_t;

#include "interpolate.hxx"
#include "inc_push.c"
#include "inc_curr.c"
#include "inc_step.c"

// ======================================================================

template<typename C>
struct PushParticles1vb
{
  using Mparticles = typename C::Mparticles;
  using Mfields = typename C::Mfields;
  using mparticles_t = PscMparticles<Mparticles>;
  using PushOne_t = PushOne;
  
  static void push_mprts(Mparticles& mprts, Mfields& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      auto& prts = mprts[p];

      flds.zero(JXI, JXI + 3);

      typename C::curr_cache_t curr_cache(flds);

      unsigned int n_prts = prts.size();
      for (int n = 0; n < n_prts; n++) {
	PushOne_t::push<C>(prts, n, flds, curr_cache);
      }
    }
  }

  static void stagger_mprts(Mparticles& mprts, Mfields& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      auto& prts = mprts[p];

      flds.zero(JXI, JXI + 3);

      unsigned int n_prts = prts.size();
      for (int n = 0; n < n_prts; n++) {
	PushOne::stagger<C>(prts, n, flds);
      }
    }
  }
};

template<typename C>
void push_p_ops<C>::push_mprts(typename C::Mparticles& mprts,
			       typename C::Mfields& mflds)
{
  PushParticles1vb<C>::push_mprts(mprts, mflds);
}

template<typename C>
void push_p_ops<C>::push_mprts(struct psc_mparticles *mprts,
			       struct psc_mfields *mflds_base)
{
  auto mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  auto mp = mparticles_t(mprts);
  PushParticles1vb<C>::push_mprts(*mp.sub(), *mf.sub());
  mf.put_as(mflds_base, JXI, JXI+3);
}

template<typename C>
void push_p_ops<C>::stagger_mprts(struct psc_mparticles *mprts,
				  struct psc_mfields *mflds_base)
{
  auto mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  auto mp = mparticles_t(mprts);
  PushParticles1vb<C>::stagger_mprts(*mp.sub(), *mf.sub());
  mf.put_as(mflds_base, JXI, JXI+3);
}

template struct push_p_ops<push_p_conf>;
