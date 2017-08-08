
#define PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end)	      \
  for (particle_iter_t prt_iter = prt_begin;			      \
       !particle_iter_equal(prt_iter, prt_end);			      \
       prt_iter = particle_iter_next(prt_iter))

static inline void
particle_range_resize(particle_range_t *prts, unsigned int n)
{
  struct psc_mparticles *mprts = (struct psc_mparticles *) prts->end.mprts;
  psc_mparticles_resize_patch(mprts, prts->end.p, n);
  prts->end.n = n;
}
