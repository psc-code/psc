
// ======================================================================
// particle_iter_t

typedef struct {
  int n;
  int p;
  const struct psc_mparticles *mprts;
} particle_iter_t;

static inline bool
particle_iter_equal(particle_iter_t iter, particle_iter_t iter2)
{
  assert(iter.mprts == iter2.mprts);
  return iter.n == iter2.n && iter.p == iter2.p;
}

static inline particle_iter_t
particle_iter_next(particle_iter_t iter)
{
  return (particle_iter_t) {
    .n    = iter.n + 1,
    .p    = iter.p,
    .mprts = iter.mprts,
  };
}

static inline particle_t *
particle_iter_deref(particle_iter_t iter)
{
  // FIXME, shouldn't have to cast away const
  return mparticles_get_one((struct psc_mparticles *) iter.mprts, iter.p, iter.n);
}

static inline particle_t *
particle_iter_at(particle_iter_t iter, int m)
{
  // FIXME, shouldn't have to cast away const
  return mparticles_get_one((struct psc_mparticles *) iter.mprts, iter.p, iter.n + m);
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
  return (particle_range_t) {
    .begin = { .n = 0                                       , .p = p, .mprts = mprts },
    .end   = { .n = psc_mparticles_n_prts_by_patch(mprts, p), .p = p, .mprts = mprts },
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
  struct psc_mparticles *mprts = (struct psc_mparticles *) prts->end.mprts;
  psc_mparticles_resize_patch(mprts, prts->end.p, n);
  prts->end.n = n;
}
