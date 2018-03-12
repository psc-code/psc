
#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_mparticles_double_ ## x
#define psc_mparticles_sub psc_mparticles_double
using mparticles_t = PscMparticlesDouble;

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_mparticles_single_ ## x
#define psc_mparticles_sub psc_mparticles_single
using mparticles_t = PscMparticlesSingle;

#elif PSC_PARTICLES_AS_FORTRAN

#define PFX(x) psc_mparticles_fortran_ ## x
#define psc_mparticles_sub psc_mparticles_fortran
using mparticles_t = PscMparticlesFortran;

#endif

template<typename MP_FROM, typename MP_TO>
struct Convert
{
  using MparticlesFrom = MP_FROM;
  using MparticlesTo = MP_TO;
  using particle_from_t = typename MP_FROM::particle_t;
  using particle_to_t = typename MP_TO::particle_t;

  particle_to_t operator()(const particle_from_t& prt_from)
  {
    particle_to_t prt_to;
    
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

template<typename Convert>
void psc_mparticles_copy_to(struct psc_mparticles *mprts_from_,
			    struct psc_mparticles *mprts_to_, unsigned int flags)
{
  using mparticles_to_t = PscMparticles<typename Convert::MparticlesTo>;
  using mparticles_from_t = PscMparticles<typename Convert::MparticlesFrom>;
  auto mprts_to = mparticles_to_t{mprts_to_};
  auto mprts_from = mparticles_from_t{mprts_from_};

  Convert convert;
  int n_patches = mprts_to->n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from->get_size_all(n_prts_by_patch);
  mprts_to->reserve_all(n_prts_by_patch);
  mprts_to->resize_all(n_prts_by_patch);
  
  for (int p = 0; p < n_patches; p++) {
    auto& prts_from = mprts_from[p];
    auto& prts_to = mprts_to[p];
    int n_prts = prts_from.size();
    for (int n = 0; n < n_prts; n++) {
      prts_to[n] = convert(prts_from[n]);
    }
  }
}

template<typename Convert>
void psc_mparticles_copy_from(struct psc_mparticles *mprts_to_,
			      struct psc_mparticles *mprts_from_, unsigned int flags)
{
  psc_mparticles_copy_to<Convert>(mprts_from_, mprts_to_, flags);
}


