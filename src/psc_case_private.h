
#ifndef PSC_CASE_PRIVATE_H
#define PSC_CASE_PRIVATE_H

#include <psc_case.h>

struct psc_case {
  struct mrc_obj obj;
  bool seed_by_time;

  struct psc *psc;
};

struct psc_case_ops {
  MRC_SUBCLASS_OPS(struct psc_case);
  void (*init_npt)(struct psc_case *_case, int kind, double x[3],
		   struct psc_particle_npt *npt);
  void (*init_field)(struct psc_case *_case, mfields_base_t *flds);
};

extern struct psc_case_ops psc_case_harris_ops;
extern struct psc_case_ops psc_case_test_xz_ops;
extern struct psc_case_ops psc_case_test_yz_ops;
extern struct psc_case_ops psc_case_harris_xy_ops;
extern struct psc_case_ops psc_case_langmuir_ops;
extern struct psc_case_ops psc_case_wakefield_ops;
extern struct psc_case_ops psc_case_thinfoil_ops;
extern struct psc_case_ops psc_case_foils_ops;
extern struct psc_case_ops psc_case_curvedfoil_ops;
extern struct psc_case_ops psc_case_singlepart_ops;
extern struct psc_case_ops psc_case_collisions_ops;
extern struct psc_case_ops psc_case_cone_ops;
extern struct psc_case_ops psc_case_microsphere_ops;

// ======================================================================

#define psc_case_ops(c) ((struct psc_case_ops *)((c)->obj.ops))

#endif
