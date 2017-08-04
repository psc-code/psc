
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

#include <mrc_obj.h>
#include <assert.h>

// ----------------------------------------------------------------------
// psc_particles class

MRC_CLASS_DECLARE(psc_particles, struct psc_particles);

void psc_particles_reorder(struct psc_particles *prts);

// ----------------------------------------------------------------------
// psc_mparticles class

struct psc_mparticles {
  struct mrc_obj obj;
  struct psc_particles **prts;
  int nr_patches;
  struct mrc_domain *domain;
  int *nr_particles_by_patch;
  unsigned int flags;          // flags, like MP_NEED_CELL_OFFSETS, ...
};

MRC_CLASS_DECLARE(psc_mparticles, struct psc_mparticles);

struct psc_mparticles_ops {
  MRC_SUBCLASS_OPS(struct psc_mparticles);
  void (*setup_internals)(struct psc_mparticles *mprts);
  unsigned int (*get_nr_particles)(struct psc_mparticles *mprts);
  void (*update_n_part)(struct psc_mparticles *mprts);
};

#define psc_mparticles_ops(mp) ((struct psc_mparticles_ops *) ((mp)->obj.ops))

typedef void (*psc_mparticles_copy_to_func_t)(int p, struct psc_mparticles *,
					      struct psc_mparticles *,
					      unsigned int);
typedef void (*psc_mparticles_copy_from_func_t)(int p, struct psc_mparticles *,
						struct psc_mparticles *,
						unsigned int);

#define MP_DONT_COPY (0x1)
#define MP_NEED_BLOCK_OFFSETS (0x0100)
#define MP_NEED_CELL_OFFSETS  (0x0200)
#define MP_BLOCKSIZE_MASK     (0x7000)
#define MP_BLOCKSIZE_1X1X1    (0x1000)
#define MP_BLOCKSIZE_2X2X2    (0x2000)
#define MP_BLOCKSIZE_4X4X4    (0x3000)
#define MP_BLOCKSIZE_8X8X8    (0x4000)
#define MP_NO_CHECKERBOARD    (0x10000)

static inline struct psc_particles *
psc_mparticles_get_patch(struct psc_mparticles *mp, int p)
{
  return mp->prts[p];
}

extern struct psc_mparticles_ops psc_mparticles_fortran_ops;
extern struct psc_mparticles_ops psc_mparticles_c_ops;
extern struct psc_mparticles_ops psc_mparticles_single_ops;
extern struct psc_mparticles_ops psc_mparticles_double_ops;
extern struct psc_mparticles_ops psc_mparticles_mix_ops;
extern struct psc_mparticles_ops psc_mparticles_sse2_ops;
extern struct psc_mparticles_ops psc_mparticles_cbe_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda_ops;
extern struct psc_mparticles_ops psc_mparticles_single_by_block_ops;

void psc_mparticles_set_domain_nr_particles(struct psc_mparticles *mparticles,
					    struct mrc_domain *domain,
					    int *nr_particles_by_patch);
void psc_mparticles_set_nr_particles(struct psc_mparticles *mprts, int *n_prts_by_patch);
int  psc_mparticles_nr_particles(struct psc_mparticles *mparticles);
int  psc_mparticles_n_prts_by_patch(struct psc_mparticles *mparticles, int p);
void psc_mparticles_resize_patch(struct psc_mparticles *mparticles, int p, int n_prts);
void psc_mparticles_setup_internals(struct psc_mparticles *mparticles);
void psc_mparticles_update_n_part(struct psc_mparticles *mparticles);
struct psc_mparticles *psc_mparticles_get_as(struct psc_mparticles *mparticles_base,
					     const char *type,
					     unsigned int flags);
void psc_mparticles_put_as(struct psc_mparticles *mparticles,
			   struct psc_mparticles *mparticles_base,
			   unsigned int flags);
void psc_mparticles_check(struct psc_mparticles *mparticles);

#endif
