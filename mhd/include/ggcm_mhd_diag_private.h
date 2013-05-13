
#ifndef GGCM_MHD_DIAG_PRIVATE_H
#define GGCM_MHD_DIAG_PRIVATE_H

#include "ggcm_mhd_diag.h"

#include <mrc_fld.h>

struct ggcm_mhd_diag {
  struct mrc_obj obj;
  struct ggcm_mhd *mhd;
};

struct ggcm_mhd_diag_ops {
  MRC_SUBCLASS_OPS(struct ggcm_mhd_diag);
  void (*run)(struct ggcm_mhd_diag *);
  void (*shutdown)(struct ggcm_mhd_diag *);
};

void ggcm_mhd_diag_c_write_one_field(struct mrc_io *io, struct mrc_f3 *f, int m,
				     const char *name, float scale, int outtype,
				     float plane);
void ggcm_mhd_diag_c_write_one_f3(struct mrc_io *io, struct mrc_f3 *f,
				  int outtype, float plane);

// ----------------------------------------------------------------------

struct mrc_io * ggcm_diag_lib_create_mrc_io(MPI_Comm comm, const char *outputmode,
					    int outtype, float sheet, int rank_diagsrv);
void ggcm_diag_lib_write_openggcm_attrs(struct mrc_io *io, const char *time_str);
void ggcm_diag_lib_make_time_string(char s[80], float time, double dacttime);

// FIXME, should go into libmrc in a more generic form
int parse_float_array(const char *str, float *arr, int n);

#endif
