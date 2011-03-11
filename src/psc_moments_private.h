
#ifndef PSC_MOMENTS_PRIVATE_H
#define PSC_MOMENTS_PRIVATE_H

#include <psc_moments.h>

struct psc_moments {
  struct mrc_obj obj;
};

struct psc_moments_ops {
  MRC_OBJ_OPS;
  void (*calc_densities)(struct psc_moments *moments,
			 mfields_base_t *flds, mparticles_base_t *particles,
			 mfields_base_t *res);
  void (*calc_v)(struct psc_moments *moments,
		 mfields_base_t *flds, mparticles_base_t *particles,
		 mfields_base_t *res);
  void (*calc_vv)(struct psc_moments *moments,
		  mfields_base_t *flds, mparticles_base_t *particles,
		  mfields_base_t *res);
};

// ======================================================================

extern struct psc_moments_ops psc_moments_c_ops;
extern struct psc_moments_ops psc_moments_fortran_ops;

#define to_psc_moments(o) (container_of(o, struct psc_moments, obj))
#define psc_moments_ops(moments) ((struct psc_moments_ops *)((moments)->obj.ops))

#endif
