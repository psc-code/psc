
static inline void
psc_mparticles_copy_from(int p, struct psc_mparticles *mprts,
			 struct psc_mparticles *mprts_from, unsigned int flags,
			 void (*get_particle)(particle_t *prt, int n, struct psc_particles *prts))
{
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles *prts_from = psc_mparticles_get_patch(mprts_from, p);
  int n_prts = psc_mparticles_n_prts_by_patch(mprts_from, p);
  psc_mparticles_resize_patch(mprts, p, n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_t *prt = particles_get_one(prts, n);
    get_particle(prt, n, prts_from);
  }
}

static inline void
psc_mparticles_copy_to(int p, struct psc_mparticles *mprts,
		       struct psc_mparticles *mprts_to, unsigned int flags,
		      void (*put_particle)(particle_t *prt, int n, struct psc_particles *prts))
{
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  struct psc_particles *prts_to = psc_mparticles_get_patch(mprts_to, p);
  int n_prts = psc_mparticles_n_prts_by_patch(mprts, p);
  psc_mparticles_resize_patch(mprts_to, p, n_prts);
  for (int n = 0; n < n_prts; n++) {
    particle_t *prt = particles_get_one(prts, n);
    put_particle(prt, n, prts_to);
  }
}


