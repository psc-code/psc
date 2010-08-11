
#ifndef OUTPUT_FIELDS_H
#define OUTPUT_FIELDS_H

#include "psc.h"

enum {
  X_NE, X_NI, X_NN,
  X_JXI, X_JYI, X_JZI,
  X_EX , X_EY , X_EZ ,
  X_HX , X_HY , X_HZ ,
  X_JXEX, X_JYEY, X_JZEZ,
  X_POYX, X_POYY, X_POYZ,
  X_E2X, X_E2Y, X_E2Z,
  X_B2X, X_B2Y, X_B2Z,
  NR_EXTRA_FIELDS,
};

#define MAX_FIELDS_LIST NR_EXTRA_FIELDS

struct psc_field {
  float *data;
  int ilo[3], ihi[3];
  unsigned int size;
  const char *name;
};

struct psc_fields_list {
  int nr_flds;
  struct psc_field flds[MAX_FIELDS_LIST];
  bool *dowrite_fd; // FIXME, obsolete -- don't use
};

struct psc_output_format_ops {
  const char *name;
  const char *ext;
  void (*create)(void);
  void (*destroy)(void);
  void (*open)(struct psc_fields_list *flds, const char *filename, void **pctx);
  void (*close)(void *ctx);
  void (*write_field)(void *ctx, struct psc_field *fld);
};

extern struct psc_output_format_ops psc_output_format_ops_binary;
extern struct psc_output_format_ops psc_output_format_ops_hdf5;
extern struct psc_output_format_ops psc_output_format_ops_vtk;
extern struct psc_output_format_ops psc_output_format_ops_vtk_points;
extern struct psc_output_format_ops psc_output_format_ops_vtk_cells;

#endif
