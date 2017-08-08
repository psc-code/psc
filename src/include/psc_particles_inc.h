
static inline void
psc_mparticles_copy_from(int p, struct psc_mparticles *mprts,
			 struct psc_mparticles *mprts_from, unsigned int flags,
			 void (*get_particle)(particle_t *prt, int n, struct psc_mparticles *mprts, int p))
{
  int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);

  for (int n = 0; n < n_prts; n++) {
    particle_t *prt = mparticles_get_one(mprts, p, n);
    get_particle(prt, n, mprts_from, p);
  }
}

static inline void
psc_mparticles_copy_to(int p, struct psc_mparticles *mprts,
		       struct psc_mparticles *mprts_to, unsigned int flags,
		       void (*put_particle)(particle_t *prt, int n, struct psc_mparticles *mprts, int p))
{
  int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);

  for (int n = 0; n < n_prts; n++) {
    particle_t *prt = mparticles_get_one(mprts, p, n);
    put_particle(prt, n, mprts_to, p);
  }
}


