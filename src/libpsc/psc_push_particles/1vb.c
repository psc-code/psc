
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "1vb/psc_push_particles_1vb.h"

#include "fields.hxx"

#include "inc_params.c"
#include "inc_cache.c"
#include "inc_interpolate.c"
#include "inc_push.c"
#include "inc_curr.c"
#include "inc_step.c"

// ======================================================================

template<typename C>
static void
do_push_part_1vb_yz(fields_t flds, struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = mparticles_t(mprts)[p].range();
  unsigned int n_prts = prts.size();
  
  for (int n = 0; n < n_prts; n++) {
    push_one<C>(prts, n, flds, flds);
  }
}

template<typename C>
static void
do_stagger_part_1vb_yz(fields_t flds, struct psc_mparticles *mprts, int p)
{
  particle_range_t prts = mparticles_t(mprts)[p].range();
  unsigned int n_prts = prts.size();
  
  for (int n = 0; n < n_prts; n++) {
    stagger_one<C>(prts, n, flds);
  }
}

template<typename C>
void push_p_ops<C>::push_mprts(struct psc_push_particles *push,
			       struct psc_mparticles *mprts,
			       struct psc_mfields *mflds_base)
{
  using mfields_t = typename C::mfields_t;
  using fields_t = typename mfields_t::fields_t;
  
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  c_prm_set(ppsc);
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = mf[p];

    flds.zero(JXI, JXI + 3);
    ext_prepare_sort_before(mprts, p);
    do_push_part_1vb_yz<C>(flds, mprts, p);
  }
  mf.put_as(mflds_base, JXI, JXI+3);
}

template<typename C>
void push_p_ops<C>::stagger_mprts(struct psc_push_particles *push,
				  struct psc_mparticles *mprts,
				  struct psc_mfields *mflds_base)
{
  using mfields_t = typename C::mfields_t;
  using fields_t = typename mfields_t::fields_t;
  
  mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  c_prm_set(ppsc);
  params_1vb_set(ppsc, NULL, NULL);
  for (int p = 0; p < mprts->nr_patches; p++) {
    fields_t flds = mf[p];
    
    flds.zero(JXI, JXI + 3);
    ext_prepare_sort_before(mprts, p);
    do_stagger_part_1vb_yz<C>(flds, mprts, p);
  }
  mf.put_as(mflds_base, JXI, JXI+3);
}

using push_p_conf = push_p_config<mfields_t, opt_dim>;

template struct push_p_ops<push_p_conf>;
