
#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_mparticles_double_ ## x
#define psc_mparticles_sub psc_mparticles_double
using mparticles_t = mparticles_double_t;

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_mparticles_single_ ## x
#define psc_mparticles_sub psc_mparticles_single
using mparticles_t = mparticles_single_t;

#elif PSC_PARTICLES_AS_FORTRAN

#define PFX(x) psc_mparticles_fortran_ ## x
#define psc_mparticles_sub psc_mparticles_fortran
using mparticles_t = mparticles_fortran_t;

#endif

template<typename F>
void psc_mparticles_copy_from(struct psc_mparticles *mprts_to,
			      struct psc_mparticles *mprts_from, unsigned int flags,
			      F convert_from)
{
  for (int p = 0; p < mprts_to->nr_patches; p++) {
    mparticles_t::patch_t& prts = mparticles_t(mprts_to)[p];
    int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      convert_from(&prts[n], n, mprts_from, p);
    }
  }
}

template<typename F>
void psc_mparticles_copy_to(struct psc_mparticles *mprts_from,
			    struct psc_mparticles *mprts_to, unsigned int flags,
			    F convert_to)
{
  for (int p = 0; p < mprts_from->nr_patches; p++) {
    mparticles_t::patch_t& prts = mparticles_t(mprts_from)[p];
    int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      convert_to(&prts[n], n, mprts_to, p);
    }
  }
}


