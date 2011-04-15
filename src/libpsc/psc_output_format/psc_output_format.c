
#include "psc_output_format_private.h"

// ======================================================================
// forward to subclass

void
psc_output_format_write_fields(struct psc_output_format *format,
			       struct psc_output_fields_c *out,
			       struct psc_fields_list *list, const char *pfx)
{
  struct psc_output_format_ops *ops = psc_output_format_ops(format);
  assert(ops->write_fields);
  ops->write_fields(format, out, list, pfx);
}

// ======================================================================
// psc_output_format_init

static void
psc_output_format_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_mrc_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_binary_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_vtk_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_vtk_points_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_vtk_cells_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_vtk_binary_ops);
#ifdef HAVE_HDF5
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_hdf5_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_format, &psc_output_format_xdmf_ops);
#endif
}

// ======================================================================
// psc_output_format class

struct mrc_class_psc_output_format mrc_class_psc_output_format = {
  .name             = "psc_output_format",
  .size             = sizeof(struct psc_output_format),
  .init             = psc_output_format_init,
};

