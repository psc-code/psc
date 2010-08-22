
#ifndef OUTPUT_FIELDS_H
#define OUTPUT_FIELDS_H

#include "psc.h"

enum {
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

struct psc_fields_list {
  int nr_flds;
  fields_base_t flds[MAX_FIELDS_LIST];
  bool *dowrite_fd; // FIXME, obsolete -- don't use
};

struct psc_output_c;

struct psc_output_format_ops {
  const char *name;
  const char *ext;
  void (*create)(void);
  void (*destroy)(void);
  void (*open)(struct psc_output_c *out, struct psc_fields_list *flds,
	       const char *prefix, void **pctx);
  void (*close)(void *ctx);
  void (*write_field)(void *ctx, fields_base_t *fld);
};

extern struct psc_output_format_ops psc_output_format_ops_binary;
extern struct psc_output_format_ops psc_output_format_ops_hdf5;
extern struct psc_output_format_ops psc_output_format_ops_xdmf;
extern struct psc_output_format_ops psc_output_format_ops_vtk;
extern struct psc_output_format_ops psc_output_format_ops_vtk_points;
extern struct psc_output_format_ops psc_output_format_ops_vtk_cells;

struct psc_output_c {
  char *data_dir;
  char *output_format;
  bool output_combine;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_first, tfield_first;
  int pfield_step, tfield_step;
  bool dowrite_fd[NR_EXTRA_FIELDS];

  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  struct psc_fields_list pfd, tfd;

  struct psc_output_format_ops *format_ops;
};

char *psc_output_c_filename(struct psc_output_c *out, const char *pfx);

#endif
