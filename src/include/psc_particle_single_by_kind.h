
#ifndef PSC_PARTICLE_SINGLE_BY_KIND_H
#define PSC_PARTICLE_SINGLE_BY_KIND_H

typedef float particle_single_by_kind_real_t;
typedef struct {
  particle_single_by_kind_real_t dx[3];
  int i;
  particle_single_by_kind_real_t ux[3];
  particle_single_by_kind_real_t w;
  int kind;
} particle_single_by_kind_t;

#define psc_mparticles_single_by_kind(mprts)({				\
      assert((struct psc_mparticles_ops *) mprts->obj.ops == &psc_mparticles_single_by_kind_ops); \
      mrc_to_subobj(mprts, struct psc_mparticles_single_by_kind);	\
})

#endif
