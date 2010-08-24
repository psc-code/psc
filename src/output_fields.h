
#ifndef OUTPUT_FIELDS_H
#define OUTPUT_FIELDS_H

#include "psc.h"

#define MAX_FIELDS_LIST 30

struct psc_fields_list {
  int nr_flds;
  fields_base_t flds[MAX_FIELDS_LIST];
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
  void (*write_fields)(void *ctx, struct psc_fields_list *flds);
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
  char *output_fields;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_first, tfield_first;
  int pfield_step, tfield_step;

  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  struct psc_fields_list pfd, tfd;
  struct output_field *out_flds[MAX_FIELDS_LIST];

  struct psc_output_format_ops *format_ops;
};

void write_fields_combine(struct psc_fields_list *list, 
			  void (*write_field)(void *ctx, fields_base_t *fld),
			  void *ctx);

#endif
