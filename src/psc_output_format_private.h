
#ifndef PSC_OUTPUT_FORMAT_PRIVATE_H
#define PSC_OUTPUT_FORMAT_PRIVATE_H

#include <psc_output_format.h>

struct psc_output_format {
  struct mrc_obj obj;
};

struct psc_output_format_ops {
  MRC_SUBCLASS_OPS(struct psc_output_format);
  void (*write_fields)(struct psc_output_format *format,
		       struct psc_output_fields_c *out,
		       struct psc_fields_list *list, const char *pfx);
};

// ======================================================================

extern struct psc_output_format_ops psc_output_format_mrc_ops;
extern struct psc_output_format_ops psc_output_format_binary_ops;
extern struct psc_output_format_ops psc_output_format_vtk_ops;
extern struct psc_output_format_ops psc_output_format_vtk_points_ops;
extern struct psc_output_format_ops psc_output_format_vtk_cells_ops;
extern struct psc_output_format_ops psc_output_format_vtk_binary_ops;
extern struct psc_output_format_ops psc_output_format_hdf5_ops;
extern struct psc_output_format_ops psc_output_format_xdmf_ops;

#define psc_output_format_ops(format) ((struct psc_output_format_ops *)((format)->obj.ops))

#endif
