
#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_mparticles_double_ ## x
#define psc_mparticles_sub psc_mparticles_double
#define PARTICLE_BUF(x) psc_particle_double_buf_ ## x
using mparticles_t = mparticles_double_t;

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_mparticles_single_ ## x
#define psc_mparticles_sub psc_mparticles_single
#define PARTICLE_BUF(x) psc_particle_single_buf_ ## x
using mparticles_t = mparticles_single_t;

#elif PSC_PARTICLES_AS_FORTRAN

#define PFX(x) psc_mparticles_fortran_ ## x
#define psc_mparticles_sub psc_mparticles_fortran
#define PARTICLE_BUF(x) psc_particle_fortran_buf_ ## x
using mparticles_t = mparticles_fortran_t;

#endif

static inline void
psc_mparticles_copy_from(struct psc_mparticles *mprts,
			 struct psc_mparticles *mprts_from, unsigned int flags,
			 void (*get_particle)(particle_t *prt, int n, struct psc_mparticles *mprts, int p))
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = mparticles_t(mprts)[p].range();
    int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      get_particle(&prts[n], n, mprts_from, p);
    }
  }
}

static inline void
psc_mparticles_copy_to(struct psc_mparticles *mprts,
		       struct psc_mparticles *mprts_to, unsigned int flags,
		       void (*put_particle)(particle_t *prt, int n, struct psc_mparticles *mprts, int p))
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = mparticles_t(mprts)[p].range();
    int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      put_particle(&prts[n], n, mprts_to, p);
    }
  }
}


