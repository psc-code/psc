
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
