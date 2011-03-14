
#ifndef PSC_CASE_PRIVATE_H
#define PSC_CASE_PRIVATE_H

#include <psc_case.h>

struct _psc_case {
  struct mrc_obj obj;
  struct psc_case *Case;
};

struct _psc_case_ops {
  MRC_SUBCLASS_OPS(struct _psc_case);
  void (*init_npt)(struct _psc_case *_case, int kind, double x[3],
		   struct psc_particle_npt *npt);
  void (*init_field)(struct _psc_case *_case, mfields_base_t *flds);
};

extern struct _psc_case_ops _psc_case_harris_ops;
extern struct _psc_case_ops _psc_case_test_xz_ops;
extern struct _psc_case_ops _psc_case_test_yz_ops;

// ======================================================================

#define _psc_case_ops(c) ((struct _psc_case_ops *)((c)->obj.ops))

#endif
