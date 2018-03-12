
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

template<typename MP, typename F>
void psc_mparticles_copy_from(mparticles_t mprts_to, MP mprts_from, unsigned int flags,
			      F convert_from)
{
  int n_patches = mprts_to->n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from->get_size_all(n_prts_by_patch);
  mprts_to->reserve_all(n_prts_by_patch);
  mprts_to->resize_all(n_prts_by_patch);
  
  for (int p = 0; p < n_patches; p++) {
    mparticles_t::patch_t& prts = mprts_to[p];
    int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      prts[n] = convert_from(mprts_from, p, n);
    }
  }
}

template<typename MP, typename F>
void psc_mparticles_copy_to(mparticles_t mprts_from, MP mprts_to, unsigned int flags,
			    F convert_to)
{
  int n_patches = mprts_to->n_patches();
  uint n_prts_by_patch[n_patches];
  mprts_from->get_size_all(n_prts_by_patch);
  mprts_to->reserve_all(n_prts_by_patch);
  mprts_to->resize_all(n_prts_by_patch);
  
  for (int p = 0; p < n_patches; p++) {
    mparticles_t::patch_t& prts = mprts_from[p];
    int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      convert_to(mprts_to, p, n, prts[n]);
    }
  }
}

template<typename ConvertFrom>
void psc_mparticles_copy_to_(struct psc_mparticles *mprts,
			     struct psc_mparticles *mprts_to, unsigned int flags)
{
  using mparticles_to_t = PscMparticles<typename ConvertFrom::MparticlesTo>;
  using mparticles_from_t = PscMparticles<typename ConvertFrom::MparticlesFrom>;

  ConvertFrom convert_from;
  psc_mparticles_copy_to(mparticles_from_t{mprts}, mparticles_to_t{mprts_to},
			 flags, convert_from);
}

template<typename ConvertTo>
void psc_mparticles_copy_from_(struct psc_mparticles *mprts,
			       struct psc_mparticles *mprts_from, unsigned int flags)
{
  using mparticles_to_t = PscMparticles<typename ConvertTo::MparticlesTo>;
  using mparticles_from_t = PscMparticles<typename ConvertTo::MparticlesFrom>;

  ConvertTo convert_to;
  psc_mparticles_copy_from(mparticles_to_t{mprts}, mparticles_from_t{mprts_from},
			   flags, convert_to);
}


