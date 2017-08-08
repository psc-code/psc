
#if PSC_PARTICLES_AS_DOUBLE

#define PFX(x) psc_particles_double_ ## x
#define MPFX(x) psc_mparticles_double_ ## x
#define psc_mparticles_sub psc_mparticles_double

#elif PSC_PARTICLES_AS_SINGLE

#define PFX(x) psc_particles_single_ ## x
#define MPFX(x) psc_mparticles_single_ ## x
#define psc_mparticles_sub psc_mparticles_single

#elif PSC_PARTICLES_AS_SINGLE_BY_BLOCK

#define PFX(x) psc_particles_single_by_block_ ## x
#define MPFX(x) psc_mparticles_single_by_block_ ## x
#define psc_mparticles_sub psc_mparticles_single_by_block

#elif PSC_PARTICLES_AS_C

#define PFX(x) psc_particles_c_ ## x
#define MPFX(x) psc_mparticles_c_ ## x
#define psc_mparticles_sub psc_mparticles_c

#elif PSC_PARTICLES_AS_FORTRAN

#define PFX(x) psc_particles_fortran_ ## x
#define MPFX(x) psc_mparticles_fortran_ ## x
#define psc_mparticles_sub psc_mparticles_fortran

#endif

static inline void
psc_mparticles_copy_from(struct psc_mparticles *mprts,
			 struct psc_mparticles *mprts_from, unsigned int flags,
			 void (*get_particle)(particle_t *prt, int n, struct psc_mparticles *mprts, int p))
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);

    for (int n = 0; n < n_prts; n++) {
      particle_t *prt = mparticles_get_one(mprts, p, n);
      get_particle(prt, n, mprts_from, p);
    }
  }
}

static inline void
psc_mparticles_copy_to(struct psc_mparticles *mprts,
		       struct psc_mparticles *mprts_to, unsigned int flags,
		       void (*put_particle)(particle_t *prt, int n, struct psc_mparticles *mprts, int p))
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);

    for (int n = 0; n < n_prts; n++) {
      particle_t *prt = mparticles_get_one(mprts, p, n);
      put_particle(prt, n, mprts_to, p);
    }
  }
}


