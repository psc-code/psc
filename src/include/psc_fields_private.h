
#ifndef PSC_FIELDS_PRIVATE_H
#define PSC_FIELDS_PRIVATE_H

#include "psc_fields.h"

struct psc_fields {
  struct mrc_obj obj;
  void *data;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  int first_comp; // first component
  int p; // patch nr
  struct psc_mfields *mflds;
};

struct psc_fields_ops {
  MRC_SUBCLASS_OPS(struct psc_fields);
  void (*zero_comp)(struct psc_fields *pf, int m);
  void (*set_comp)(struct psc_fields *pf, int m, double alpha);
  void (*scale_comp)(struct psc_fields *pf, int m, double alpha);
  void (*copy_comp)(struct psc_fields *to, int mto,
		    struct psc_fields *from, int mfrom);
  void (*axpy_comp)(struct psc_fields *yf, int ym, double alpha,
		    struct psc_fields *xf, int xm);
};

#define psc_fields_ops(prts) ((struct psc_fields_ops *) ((prts)->obj.ops))

typedef void (*psc_fields_copy_to_func_t)(struct psc_fields *,
					  struct psc_fields *,
					  int, int);
typedef void (*psc_fields_copy_from_func_t)(struct psc_fields *,
					    struct psc_fields *,
					    int, int);

// ======================================================================

extern struct psc_fields_ops psc_fields_c_ops;
extern struct psc_fields_ops psc_fields_single_ops;
extern struct psc_fields_ops psc_fields_fortran_ops;
extern struct psc_fields_ops psc_fields_cuda_ops;
extern struct psc_fields_ops psc_fields_cuda2_ops;
extern struct psc_fields_ops psc_fields_acc_ops;

extern struct psc_mfields_ops psc_mfields_cuda2_ops;
extern struct psc_mfields_ops psc_mfields_acc_ops;

#endif
