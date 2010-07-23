
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

struct psc_extra_fields {
  unsigned int size;
  unsigned int naccum;
  float *all[NR_EXTRA_FIELDS];
};


void init_output_fields(struct psc_extra_fields *f);
void free_output_fields(struct psc_extra_fields *f);
void calculate_pfields(struct psc_extra_fields *p);
void accumulate_tfields(struct psc_extra_fields *p, struct psc_extra_fields *t);
void reset_fields(struct psc_extra_fields *f);

// convert accumulated values to correct temporal mean
// (divide by naccum)
void mean_tfields(struct psc_extra_fields *f);

struct psc_output_format_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*write_fields)(struct psc_extra_fields *flds, bool *dowrite_fd,
		       const char *prefix);
};

extern struct psc_output_format_ops psc_output_format_ops_binary;

#endif
