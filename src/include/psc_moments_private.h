
#ifndef PSC_MOMENTS_PRIVATE_H
#define PSC_MOMENTS_PRIVATE_H

#include <psc_moments.h>

struct psc_moments {
  struct mrc_obj obj;
  struct psc_bnd *bnd; // used for add_ghosts, etc.
};

struct psc_moments_ops {
  MRC_SUBCLASS_OPS(struct psc_moments);
  void (*calc_densities)(struct psc_moments *moments,
			 mfields_base_t *flds, mparticles_base_t *particles,
			 mfields_c_t *res);
  void (*calc_v)(struct psc_moments *moments,
		 mfields_base_t *flds, mparticles_base_t *particles,
		 mfields_c_t *res);
  void (*calc_vv)(struct psc_moments *moments,
		  mfields_base_t *flds, mparticles_base_t *particles,
		  mfields_c_t *res);
  void (*calc_photon_n)(struct psc_moments *moments,
			mphotons_t *photons, mfields_c_t *res);
};

// ======================================================================

extern struct psc_moments_ops psc_moments_c_ops;
extern struct psc_moments_ops psc_moments_1st_ops;

#define psc_moments_ops(moments) ((struct psc_moments_ops *)((moments)->obj.ops))

#endif
