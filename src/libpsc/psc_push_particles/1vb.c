
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "1vb/psc_push_particles_1vb.h"

#include "fields.hxx"

using real_t = mparticles_t::real_t;

#include "inc_params.c"
#include "inc_cache.c"
#include "interpolate.hxx"
#include "inc_push.c"
#include "inc_curr.c"
#include "inc_step.c"

// ======================================================================

template<typename C>
struct PushParticles1vb
{
  using mparticles_t = typename C::mparticles_t;
  using mfields_t = typename C::mfields_t;
  using Mparticles = typename C::mparticles_t::sub_t;
  using Mfields = typename C::mfields_t::sub_t;
  
  static void push_mprts(Mparticles& mprts, Mfields& mflds)
  {
    c_prm_set(mprts.grid());
    params_1vb_set(mprts.grid());
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      auto& prts = mprts[p];

      flds.zero(JXI, JXI + 3);
      ExtPrepareSort<mparticles_t, typename C::ext>::before(prts);

      unsigned int n_prts = prts.size();
      for (int n = 0; n < n_prts; n++) {
	push_one<C>(prts, n, flds, flds);
      }
    }
  }

  static void stagger_mprts(Mparticles& mprts, Mfields& mflds)
  {
    c_prm_set(mprts.grid());
    params_1vb_set(mprts.grid());
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      auto& prts = mprts[p];

      flds.zero(JXI, JXI + 3);
      ExtPrepareSort<mparticles_t, typename C::ext>::before(prts);

      unsigned int n_prts = prts.size();
      for (int n = 0; n < n_prts; n++) {
	stagger_one<C>(prts, n, flds);
      }
    }
  }
};

template<typename C>
void push_p_ops<C>::push_mprts(typename C::mparticles_t::sub_t& mprts,
			       typename C::mfields_t::sub_t& mflds)
{
  PushParticles1vb<C>::push_mprts(mprts, mflds);
}

template<typename C>
void push_p_ops<C>::push_mprts(struct psc_push_particles *push,
			       struct psc_mparticles *mprts,
			       struct psc_mfields *mflds_base)
{
  using mfields_t = typename C::mfields_t;
  using mparticles_t = typename C::mparticles_t;
  
  auto mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  auto mp = mparticles_t(mprts);
  PushParticles1vb<C>::push_mprts(*mp.sub(), *mf.sub());
  mf.put_as(mflds_base, JXI, JXI+3);
}

template<typename C>
void push_p_ops<C>::stagger_mprts(struct psc_push_particles *push,
				  struct psc_mparticles *mprts,
				  struct psc_mfields *mflds_base)
{
  using mfields_t = typename C::mfields_t;
  using mparticles_t = typename C::mparticles_t;
  
  auto mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  auto mp = mparticles_t(mprts);
  PushParticles1vb<C>::stagger_mprts(*mp.sub(), *mf.sub());
  mf.put_as(mflds_base, JXI, JXI+3);
}

template struct push_p_ops<push_p_conf>;
