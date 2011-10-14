
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

// ----------------------------------------------------------------------
// mparticles type

struct psc_mparticles {
  struct mrc_obj obj;
  void *data;
  int nr_patches;
  struct mrc_domain *domain;
  int *nr_particles_by_patch;
};

MRC_CLASS_DECLARE(psc_mparticles, struct psc_mparticles);

struct psc_mparticles_ops {
  MRC_SUBCLASS_OPS(struct psc_mparticles);
  int  (*nr_particles_by_patch)(struct psc_mparticles *mparticles, int p);
  void (*copy_to_c)(struct psc_mparticles *particles_base, struct psc_mparticles *particles,
		    unsigned int flags);
  void (*copy_to_fortran)(struct psc_mparticles *particles_base, struct psc_mparticles *particles,
			  unsigned int flags);
  void (*copy_to_cuda)(struct psc_mparticles *particles_base, struct psc_mparticles *particles,
		       unsigned int flags);
  void (*copy_from_c)(struct psc_mparticles *particles_base, struct psc_mparticles *particles,
		      unsigned int flags);
  void (*copy_from_fortran)(struct psc_mparticles *particles_base, struct psc_mparticles *particles,
			    unsigned int flags);
  void (*copy_from_cuda)(struct psc_mparticles *particles_base, struct psc_mparticles *particles,
			 unsigned int flags);
};

#define MP_DONT_COPY (0x1)
#define MP_NEED_BLOCK_OFFSETS (0x100)
#define MP_NEED_CELL_OFFSETS (0x200)

typedef struct psc_mparticles mparticles_c_t;
typedef struct psc_mparticles mparticles_fortran_t;
typedef struct psc_mparticles mparticles_sse2_t;
typedef struct psc_mparticles mparticles_cbe_t;
typedef struct psc_mparticles mparticles_cuda_t;
typedef struct psc_mparticles mparticles_base_t;

#include "psc_particles_fortran.h"
#include "psc_particles_c.h"
#ifdef USE_SSE2
#include "psc_particles_sse2.h"
#endif

#include "psc_particles_cbe.h"

#ifdef USE_CUDA
#include "psc_particles_cuda.h"
#endif

extern struct psc_mparticles_ops psc_mparticles_c_ops;
extern struct psc_mparticles_ops psc_mparticles_fortran_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda_ops;
extern struct psc_mparticles_ops psc_mparticles_cbe_ops;

#define psc_mparticles_ops(mp) (struct psc_mparticles_ops *) ((mp)->obj.ops)

static inline particles_c_t *
psc_mparticles_get_patch_c(mparticles_c_t *mp, int p)
{
  assert(psc_mparticles_ops(mp) == &psc_mparticles_c_ops);
  return ((particles_c_t *)mp->data) + p;
}

static inline particles_fortran_t *
psc_mparticles_get_patch_fortran(mparticles_fortran_t *mp, int p)
{
  assert(psc_mparticles_ops(mp) == &psc_mparticles_fortran_ops);
  return ((particles_fortran_t *)mp->data) + p;
}

#ifdef USE_CUDA
static inline particles_cuda_t *
psc_mparticles_get_patch_cuda(mparticles_cuda_t *mp, int p)
{
  assert(psc_mparticles_ops(mp) == &psc_mparticles_cuda_ops);
  return ((particles_cuda_t *)mp->data) + p;
}
#endif

static inline particles_cbe_t *
psc_mparticles_get_patch_cbe(mparticles_cbe_t *mp, int p)
{
  assert(psc_mparticles_ops(mp) == &psc_mparticles_cbe_ops);
  return ((particles_cbe_t *)mp->data) + p;
}

void psc_mparticles_set_domain_nr_particles(struct psc_mparticles *mparticles,
					    struct mrc_domain *domain,
					    int *nr_particles_by_patch);
int  psc_mparticles_nr_particles_by_patch(struct psc_mparticles *mparticles, int p);
mparticles_c_t *psc_mparticles_get_c(struct psc_mparticles *mparticles_base,
				     unsigned int flags);
void psc_mparticles_put_c(mparticles_c_t *mparticles, struct psc_mparticles *particles_base);
mparticles_fortran_t *psc_mparticles_get_fortran(struct psc_mparticles *mparticles_base,
						 unsigned int flags);
void psc_mparticles_put_fortran(mparticles_fortran_t *mparticles, struct psc_mparticles *particles_base);
mparticles_cuda_t *psc_mparticles_get_cuda(struct psc_mparticles *mparticles_base,
					   unsigned int flags);
void psc_mparticles_put_cuda(mparticles_cuda_t *mparticles, struct psc_mparticles *particles_base);

void _psc_mparticles_fortran_copy_to_c(mparticles_fortran_t *, mparticles_c_t *, unsigned int flags);
void _psc_mparticles_fortran_copy_from_c(mparticles_fortran_t *, mparticles_c_t *, unsigned int flags);
void _psc_mparticles_cuda_copy_to_c(mparticles_cuda_t *, mparticles_c_t *, unsigned int flags);
void _psc_mparticles_cuda_copy_from_c(mparticles_cuda_t *, mparticles_c_t *, unsigned int flags);

#endif
