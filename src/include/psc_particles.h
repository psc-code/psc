
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

#include <mrc_obj.h>
#include <assert.h>

// ----------------------------------------------------------------------
// particles type

MRC_CLASS_DECLARE(psc_particles, struct psc_particles);

// ----------------------------------------------------------------------
// mparticles type

struct psc_mparticles {
  struct mrc_obj obj;
  void **patches;
  int nr_patches;
  struct mrc_domain *domain;
  int *nr_particles_by_patch;
  unsigned int flags;          // flags, like MP_NEED_CELL_OFFSETS, ...
};

MRC_CLASS_DECLARE(psc_mparticles, struct psc_mparticles);

typedef void (*psc_mparticles_copy_to_func_t)(int p, struct psc_mparticles *,
					      struct psc_mparticles *,
					      unsigned int);
typedef void (*psc_mparticles_copy_from_func_t)(int p, struct psc_mparticles *,
						struct psc_mparticles *,
						unsigned int);

struct psc_mparticles_ops {
  MRC_SUBCLASS_OPS(struct psc_mparticles);
  int  (*nr_particles_by_patch)(struct psc_mparticles *mparticles, int p);
  void *(*alloc_patch)(int p, int n_part, unsigned int flags);
  void (*free_patch)(int p, void *pp);
  int size_of_particles_t;
};

#define psc_mparticles_ops(mp) ((struct psc_mparticles_ops *) ((mp)->obj.ops))

#define MP_DONT_COPY (0x1)
#define MP_NEED_BLOCK_OFFSETS (0x0100)
#define MP_NEED_CELL_OFFSETS  (0x0200)
#define MP_BLOCKSIZE_MASK     (0x7000)
#define MP_BLOCKSIZE_1X1X1    (0x1000)
#define MP_BLOCKSIZE_2X2X2    (0x2000)
#define MP_BLOCKSIZE_4X4X4    (0x3000)
#define MP_BLOCKSIZE_8X8X8    (0x4000)
#define MP_NO_CHECKERBOARD    (0x10000)
#define MP_INTERNAL_PARTICLE_EXCHANGE (0x20000)

typedef struct psc_mparticles mparticles_base_t;

#define MAKE_MPARTICLES_TYPE(type)					\
typedef struct psc_mparticles mparticles_##type##_t;			\
									\
extern struct psc_mparticles_ops psc_mparticles_##type##_ops;		\
									\
static inline struct psc_particles_##type *				\
psc_mparticles_get_patch_##type(struct psc_mparticles *mp, int p)	\
{									\
  assert(psc_mparticles_ops(mp) == &psc_mparticles_##type##_ops);	\
  return (struct psc_particles_##type *) mp->patches[p];		\
}									\
									\
struct psc_mparticles *						        \
psc_mparticles_get_##type(struct psc_mparticles *mp_base,		\
			  unsigned int flags);				\
void psc_mparticles_put_##type(struct psc_mparticles *mp,		\
			       struct psc_mparticles *mp_base,		\
			       unsigned int flags);			\

static inline struct psc_particles *
psc_mparticles_get_patch(struct psc_mparticles *mp, int p)
{
  return mp->patches[p];
}

MAKE_MPARTICLES_TYPE(fortran)
MAKE_MPARTICLES_TYPE(c)
MAKE_MPARTICLES_TYPE(single)
MAKE_MPARTICLES_TYPE(double)
#ifdef xUSE_SSE2
MAKE_MPARTICLES_TYPE(sse2)
#endif
MAKE_MPARTICLES_TYPE(cbe)
#ifdef USE_CUDA
MAKE_MPARTICLES_TYPE(cuda)
#endif

void psc_mparticles_set_domain_nr_particles(struct psc_mparticles *mparticles,
					    struct mrc_domain *domain,
					    int *nr_particles_by_patch);
int  psc_mparticles_nr_particles(struct psc_mparticles *mparticles);
int  psc_mparticles_nr_particles_by_patch(struct psc_mparticles *mparticles, int p);
struct psc_mparticles *psc_mparticles_get_as(struct psc_mparticles *mparticles_base,
					     const char *type,
					     unsigned int flags);
void psc_mparticles_put_as(struct psc_mparticles *mparticles,
			   struct psc_mparticles *mparticles_base,
			   unsigned int flags);
void psc_mparticles_check(struct psc_mparticles *mparticles);

#endif
