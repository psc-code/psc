
#ifndef OUTPUT_FIELDS_H
#define OUTPUT_FIELDS_H

#include "psc.h"
#include "psc_output_fields_c.h"

struct psc_output_c;

struct psc_output_format_ops {
  const char *name;
  void (*create)(void);
  void (*destroy)(void);
  void (*write_fields)(struct psc_output_fields_c *out, struct psc_fields_list *flds,
		       const char *prefix);
};

extern struct psc_output_format_ops psc_output_format_ops_binary;
extern struct psc_output_format_ops psc_output_format_ops_hdf5;
extern struct psc_output_format_ops psc_output_format_ops_xdmf;
extern struct psc_output_format_ops psc_output_format_ops_vtk;
extern struct psc_output_format_ops psc_output_format_ops_vtk_points;
extern struct psc_output_format_ops psc_output_format_ops_vtk_cells;
extern struct psc_output_format_ops psc_output_format_ops_vtk_binary;
extern struct psc_output_format_ops psc_output_format_ops_mrc;

void write_fields_combine(struct psc_fields_list *list, 
			  void (*write_field)(void *ctx, fields_base_t *fld),
			  void *ctx);

#endif
