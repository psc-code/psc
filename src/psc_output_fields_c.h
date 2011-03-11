
#ifndef PSC_OUTPUT_FIELDS_C_H
#define PSC_OUTPUT_FIELDS_C_H

#include "psc_output_fields_private.h"

#define MAX_FIELDS_LIST 30

struct psc_fields_list {
  int nr_flds;
  mfields_base_t flds[MAX_FIELDS_LIST];
};

struct psc_output_fields_c {
  char *data_dir;
  char *output_format;
  char *output_fields;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_first, tfield_first;
  int pfield_step, tfield_step;
  int rn[3];
  int rx[3];
	
  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  struct psc_fields_list pfd, tfd;
  struct output_field *out_flds[MAX_FIELDS_LIST];

  struct psc_output_format_ops *format_ops;
};

#endif
