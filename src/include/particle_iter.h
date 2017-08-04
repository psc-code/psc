
// ======================================================================
// particle_iter_t

typedef struct {
  int n;
  const struct psc_particles *prts;
} particle_iter_t;

static inline bool
particle_iter_equal(particle_iter_t iter, particle_iter_t iter2)
{
  assert(iter.prts == iter2.prts);
  return iter.n == iter2.n;
}

static inline particle_iter_t
particle_iter_next(particle_iter_t iter)
{
  return (particle_iter_t) {
    .n    = iter.n + 1,
    .prts = iter.prts,
  };
}

static inline particle_t *
particle_iter_deref(particle_iter_t iter)
{
  // FIXME, shouldn't have to cast away const
  return particles_get_one((struct psc_particles *) iter.prts, iter.n);
}

static inline particle_t *
particle_iter_at(particle_iter_t iter, int m)
{
  // FIXME, shouldn't have to cast away const
  return particles_get_one((struct psc_particles *) iter.prts, iter.n + m);
}

#define PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end)	      \
  for (particle_iter_t prt_iter = prt_begin;			      \
       !particle_iter_equal(prt_iter, prt_end);			      \
       prt_iter = particle_iter_next(prt_iter))

// ----------------------------------------------------------------------
// particle_range_t

typedef struct {
  particle_iter_t begin;
  particle_iter_t end;
} particle_range_t;

static inline particle_range_t
particle_range_mprts(struct psc_mparticles *mprts, int p)
{
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
  return (particle_range_t) {
    .begin = { .n = 0                       , .prts = prts },
    .end   = { .n = psc_particles_size(prts), .prts = prts },
  };
}

static inline unsigned int
particle_range_size(particle_range_t prts)
{
  return prts.end.n - prts.begin.n;
}

static inline void
particle_range_resize(particle_range_t *prts, unsigned int n)
{
  struct psc_particles *_prts = (struct psc_particles *) prts->end.prts;
  psc_particles_resize(_prts, n);
  prts->end.n = n;
}
