
#ifndef PSC_BND_PARTICLES_PRIVATE_H
#define PSC_BND_PARTICLES_PRIVATE_H

#include <psc_bnd_particles.h>

struct psc_bnd_particles {
  struct mrc_obj obj;
  struct psc *psc;
  struct ddc_particles *ddcp;

  // for open b.c.
  double time_relax;
  bool first_time;
  struct psc_bnd *flds_bnd;
  struct psc_output_fields_item *item_nvt;
  struct psc_mfields *mflds_nvt_av;
  struct psc_mfields *mflds_nvt_last;
  struct psc_mfields *mflds_n_in;
};

struct psc_bnd_particles_ops {
  MRC_SUBCLASS_OPS(struct psc_bnd_particles);

  void (*unsetup)(struct psc_bnd_particles *bnd);
  void (*exchange_particles)(struct psc_bnd_particles *bnd, mparticles_base_t *particles);
  void (*exchange_particles_prep)(struct psc_bnd_particles *bnd, struct psc_particles *prts);
  void (*exchange_particles_post)(struct psc_bnd_particles *bnd, struct psc_particles *prts);
  void (*exchange_mprts_prep)(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts);
  void (*exchange_mprts_post)(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts);
  void (*open_calc_moments)(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts);
};

#define psc_bnd_particles_ops(bnd) ((struct psc_bnd_particles_ops *)((bnd)->obj.ops))

// ======================================================================

extern struct psc_bnd_particles_ops psc_bnd_particles_auto_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_c_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_single_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_double_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_double_omp_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_single2_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_fortran_ops;
extern struct psc_bnd_particles_ops psc_bnd_particles_cuda_ops;

#endif
