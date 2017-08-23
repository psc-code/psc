
#ifndef PSC_OUTPUT_FIELDS_C_H
#define PSC_OUTPUT_FIELDS_C_H

#include "psc_output_fields_private.h"
#include "psc_output_fields_item.h"

// ----------------------------------------------------------------------

enum {
  IO_TYPE_PFD,
  IO_TYPE_TFD,
  NR_IO_TYPES,
};

#define MAX_FIELDS_LIST 50

struct psc_fields_list {
  int nr_flds;
  struct psc_mfields *flds[MAX_FIELDS_LIST];
};

struct psc_output_fields_c {
  char *data_dir;
  char *output_fields;
  char *pfd_s;
  char *tfd_s;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_first, tfield_first;
  int pfield_step, tfield_step;
  int tfield_length;
  int tfield_every;
  int rn[3];
  int rx[3];
	
  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  struct psc_fields_list pfd, tfd;
  struct psc_output_fields_item *item[MAX_FIELDS_LIST];
  struct mrc_io *ios[NR_IO_TYPES];

  struct psc_bnd *bnd;
};

// ======================================================================

extern struct _psc_output_format_ops psc_output_format_ops_mrc;
extern struct _psc_output_format_ops psc_output_format_ops_binary;
extern struct _psc_output_format_ops psc_output_format_ops_hdf5;
extern struct _psc_output_format_ops psc_output_format_ops_xdmf;
extern struct _psc_output_format_ops psc_output_format_ops_vtk;
extern struct _psc_output_format_ops psc_output_format_ops_vtk_points;
extern struct _psc_output_format_ops psc_output_format_ops_vtk_cells;
extern struct _psc_output_format_ops psc_output_format_ops_vtk_binary;

void write_fields_combine(struct psc_fields_list *list, 
			  void (*write_field)(void *ctx, struct psc_mfields *fld),
			  void *ctx);



#endif
