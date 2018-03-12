
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

#include <mrc_obj.h>
#include <psc_bits.h>
#include <assert.h>

#include "grid.hxx"
#include "particles_traits.hxx"

// ----------------------------------------------------------------------
// psc_mparticles class

struct psc_mparticles {
  struct mrc_obj obj;
  const Grid_t* grid;
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
};

#define psc_mparticles_ops(mp) ((struct psc_mparticles_ops *) ((mp)->obj.ops))

typedef void (*psc_mparticles_copy_func_t)(struct psc_mparticles *,
					   struct psc_mparticles *,
					   unsigned int);

#define MP_DONT_COPY (0x1)
#define MP_DONT_RESIZE (0x2)
#define MP_BLOCKSIZE_MASK     (0x7000)
#define MP_BLOCKSIZE_1X1X1    (0x1000)
#define MP_BLOCKSIZE_2X2X2    (0x2000)
#define MP_BLOCKSIZE_4X4X4    (0x3000)
#define MP_BLOCKSIZE_8X8X8    (0x4000)


void psc_mparticles_check(struct psc_mparticles *mprts);

#endif
