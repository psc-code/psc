
// ======================================================================
// particle_iter_t

typedef struct {
  int n;
  const struct psc_particles *prts;
} particle_iter_t;

static inline particle_iter_t
particle_iter_begin_prts(const struct psc_particles *prts)
{
  return (particle_iter_t) {
    .n    = 0,
    .prts = prts,
  };
}

static inline particle_iter_t
particle_iter_end_prts(const struct psc_particles *prts)
{
  return (particle_iter_t) {
    .n    = prts->n_part,
    .prts = prts,
  };
}

static inline particle_iter_t
particle_iter_init(particle_iter_t iter)
{
  return (particle_iter_t) {
    .n = iter.n,
    .prts = iter.prts,
  };
}

static inline void
particle_iter_copy(particle_iter_t *iter_to, const particle_iter_t *iter_from)
{
  *iter_to = particle_iter_init(*iter_from);
}

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

#define PARTICLE_ITER_LOOP(prt_iter, prt_begin, prt_end)	      \
  for (particle_iter_t prt_iter = particle_iter_init(prt_begin);      \
       !particle_iter_equal(prt_iter, prt_end);			      \
       prt_iter = particle_iter_next(prt_iter))

