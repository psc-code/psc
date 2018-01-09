
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

#include <mrc_obj.h>
#include <assert.h>

BEGIN_C_DECLS

// ----------------------------------------------------------------------
// psc_mparticles class

struct psc_mparticles {
  struct mrc_obj obj;
  int nr_patches;
  unsigned int flags;          // flags, like MP_NEED_CELL_OFFSETS, ...
};

MRC_CLASS_DECLARE(psc_mparticles, struct psc_mparticles);

typedef struct psc_particle_inject {
  double x[3];
  double u[3];
  double w;
  int kind;
} particle_inject_t;

struct psc_mparticles_ops {
  MRC_SUBCLASS_OPS(struct psc_mparticles);
  void (*setup_internals)(struct psc_mparticles *mprts);
  unsigned int (*get_nr_particles)(struct psc_mparticles *mprts);
  void (*reserve_all)(struct psc_mparticles *mprts, int *n_prts_by_patch);
  void (*resize_all)(struct psc_mparticles *mprts, int *n_prts_by_patch);
  void (*get_size_all)(struct psc_mparticles *mprts, int *n_prts_by_patch);
  void (*inject)(struct psc_mparticles *mprts, int p,
		 const struct psc_particle_inject *prt);
};

#define psc_mparticles_ops(mp) ((struct psc_mparticles_ops *) ((mp)->obj.ops))

typedef void (*psc_mparticles_copy_func_t)(struct psc_mparticles *,
					   struct psc_mparticles *,
					   unsigned int);

#define MP_DONT_COPY (0x1)
#define MP_DONT_RESIZE (0x2)
#define MP_NEED_BLOCK_OFFSETS (0x0100)
#define MP_NEED_CELL_OFFSETS  (0x0200)
#define MP_BLOCKSIZE_MASK     (0x7000)
#define MP_BLOCKSIZE_1X1X1    (0x1000)
#define MP_BLOCKSIZE_2X2X2    (0x2000)
#define MP_BLOCKSIZE_4X4X4    (0x3000)
#define MP_BLOCKSIZE_8X8X8    (0x4000)
#define MP_NO_CHECKERBOARD    (0x10000)

extern struct psc_mparticles_ops psc_mparticles_fortran_ops;
extern struct psc_mparticles_ops psc_mparticles_single_ops;
extern struct psc_mparticles_ops psc_mparticles_double_ops;
extern struct psc_mparticles_ops psc_mparticles_sse2_ops;
extern struct psc_mparticles_ops psc_mparticles_cbe_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda_ops;
extern struct psc_mparticles_ops psc_mparticles_single_by_block_ops;
extern struct psc_mparticles_ops psc_mparticles_cuda2_ops;
extern struct psc_mparticles_ops psc_mparticles_acc_ops;
extern struct psc_mparticles_ops psc_mparticles_vpic_ops;
extern struct psc_mparticles_ops psc_mparticles_single_by_kind_ops;


int  psc_mparticles_nr_particles(struct psc_mparticles *mprts);
void psc_mparticles_get_size_all(struct psc_mparticles *mprts, int *n_prts_by_patch);
void psc_mparticles_setup_internals(struct psc_mparticles *mprts);
void psc_mparticles_reserve_all(struct psc_mparticles *mprts, int *n_prts_by_patch);
void psc_mparticles_resize_all(struct psc_mparticles *mprts, int *n_prts_by_patch);

void psc_mparticles_inject(struct psc_mparticles *mprts, int p,
			   const struct psc_particle_inject *prt);

struct psc_mparticles *psc_mparticles_get_as(struct psc_mparticles *mprts_base,
					     const char *type,
					     unsigned int flags);
void psc_mparticles_put_as(struct psc_mparticles *mprts,
			   struct psc_mparticles *mprts_base,
			   unsigned int flags);
void psc_mparticles_check(struct psc_mparticles *mprts);

END_C_DECLS

#endif
