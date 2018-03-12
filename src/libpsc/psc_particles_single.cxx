
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

template<typename MP_FROM, typename MP_TO>
struct ConvertFrom
{
  using MparticlesFrom = MP_FROM;
  using MparticlesTo = MP_TO;
  using particle_from_t = typename MP_FROM::particle_t;
  using mparticles_to_t = PscMparticles<MP_TO>;

  void operator()(mparticles_to_t mprts_to, int p, int n, const particle_t& prt_from)
  {
    auto& prt_to = mprts_to[p][n];
    
    prt_to.xi       = prt_from.xi;
    prt_to.yi       = prt_from.yi;
    prt_to.zi       = prt_from.zi;
    prt_to.pxi      = prt_from.pxi;
    prt_to.pyi      = prt_from.pyi;
    prt_to.pzi      = prt_from.pzi;
    prt_to.qni_wni_ = prt_from.qni_wni_;
    prt_to.kind_    = prt_from.kind_;
  }
};

template<typename MP_TO, typename MP_FROM>
struct ConvertTo
{
  using MparticlesFrom = MP_FROM;
  using MparticlesTo = MP_TO;
  using particle_to_t = typename MP_TO::particle_t;
  using mparticles_from_t = PscMparticles<MP_FROM>;

  particle_to_t operator()(mparticles_from_t mprts_from, int p, int n)
  {
    particle_to_t prt_to;
    const auto& prt_from = mprts_from[p][n];
    
    prt_to.xi       = prt_from.xi;
    prt_to.yi       = prt_from.yi;
    prt_to.zi       = prt_from.zi;
    prt_to.pxi      = prt_from.pxi;
    prt_to.pyi      = prt_from.pyi;
    prt_to.pzi      = prt_from.pzi;
    prt_to.qni_wni_ = prt_from.qni_wni_;
    prt_to.kind_    = prt_from.kind_;
    
    return prt_to;
  }
};

template<>
mrc_obj_method MparticlesSingle::methods[] = {
  MRC_OBJ_METHOD("copy_to_double"  , (psc_mparticles_copy_to_<ConvertFrom<MparticlesSingle, MparticlesDouble>>)),
  MRC_OBJ_METHOD("copy_from_double", (psc_mparticles_copy_from_<ConvertTo<MparticlesSingle, MparticlesDouble>>)),
  {}
};

#include "psc_particles_common.cxx"

