
#ifndef OUTPUT_FIELDS_H
#define OUTPUT_FIELDS_H

#include "psc.h"

enum {
  X_EX , X_EY , X_EZ ,
  X_BX , X_BY , X_BZ ,
  X_JXI, X_JYI, X_JZI,
  X_JXEX, X_JYEY, X_JZEZ,
  X_POYX, X_POYY, X_POYZ,
  X_E2X, X_E2Y, X_E2Z,
  X_B2X, X_B2Y, X_B2Z,
  NR_EXTRA_FIELDS,
};

#define MAX_FIELDS_LIST NR_EXTRA_FIELDS

struct psc_field {
  float *data;
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
  void (*create)(void);
  void (*destroy)(void);
  void (*write_fields)(struct psc_fields_list *flds, const char *prefix);
};

extern struct psc_output_format_ops psc_output_format_ops_binary;
extern struct psc_output_format_ops psc_output_format_ops_hdf5;

#endif
