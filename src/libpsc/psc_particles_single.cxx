
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_particles_inc.h"
#include "psc_particles_double.h"

#if 0
static void _mrc_unused // FIXME
psc_particles_single_reorder(struct psc_mparticles *mprts)
{
  struct psc_particles_single *sub = psc_particles_single(prts);

  if (!sub->need_reorder) {
    return;
  }
f
  int n_prts = patch->n_prts;
  for (int n = 0; n < n_prts; n++) {
    sub->particles_alt[n] = sub->particles[sub->b_ids[n]];
  }
  
  // swap in alt array
  particle_single_t *tmp = sub->particles;
  sub->particles = sub->particles_alt;
  sub->particles_alt = tmp;
  sub->need_reorder = false;
}
#endif

// ======================================================================
// psc_mparticles: subclass "single"

// ----------------------------------------------------------------------
// conversion to/from "double"

struct ConvertFromSingle
{
  void operator()(particle_single_t *prt, int n, struct psc_mparticles *mprts_dbl, int p)
  {
    particle_double_t *prt_dbl = &mparticles_double_t(mprts_dbl)[p][n];
    
    prt_dbl->xi      = prt->xi;
    prt_dbl->yi      = prt->yi;
    prt_dbl->zi      = prt->zi;
    prt_dbl->pxi     = prt->pxi;
    prt_dbl->pyi     = prt->pyi;
    prt_dbl->pzi     = prt->pzi;
    prt_dbl->qni_wni = prt->qni_wni;
    prt_dbl->kind_   = prt->kind_;
  }
};

struct ConvertToSingle
{
  void operator()(particle_single_t *prt, int n, struct psc_mparticles *mprts_dbl, int p)
  {
    particle_double_t *prt_dbl = &mparticles_double_t(mprts_dbl)[p][n];
    
    prt->xi      = prt_dbl->xi;
    prt->yi      = prt_dbl->yi;
    prt->zi      = prt_dbl->zi;
    prt->pxi     = prt_dbl->pxi;
    prt->pyi     = prt_dbl->pyi;
    prt->pzi     = prt_dbl->pzi;
    prt->qni_wni = prt_dbl->qni_wni;
    prt->kind_   = prt_dbl->kind_;
  }
};

static void
psc_mparticles_single_copy_to_double(struct psc_mparticles *mprts,
				     struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  ConvertFromSingle convert_from_single;
  psc_mparticles_copy_to(mparticles_t(mprts), mprts_dbl, flags, convert_from_single);
}

static void
psc_mparticles_single_copy_from_double(struct psc_mparticles *mprts,
				       struct psc_mparticles *mprts_dbl, unsigned int flags)
{
  ConvertToSingle convert_to_single;
  psc_mparticles_copy_from(mparticles_t(mprts), mprts_dbl, flags, convert_to_single);
}

static struct mrc_obj_method psc_mparticles_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_double"  , psc_mparticles_single_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", psc_mparticles_single_copy_from_double),
  {}
};

#include "psc_particles_common.cxx"

