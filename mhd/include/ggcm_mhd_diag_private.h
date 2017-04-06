
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
  void (*run_now)(struct ggcm_mhd_diag *diag, struct mrc_fld *fld,
		  int diag_type, int itdia);
  void (*shutdown)(struct ggcm_mhd_diag *);

  void (*mod_register)(struct ggcm_mhd_diag *mhd_diag, struct mrc_mod *mod);
};

void ggcm_mhd_diag_c_write_one_field(struct mrc_io *io, struct mrc_fld *f, int m,
				     const char *name, float scale, int outtype,
				     float plane);
void ggcm_mhd_diag_c_write_one_fld(struct mrc_io *io, struct mrc_fld *f,
				   int outtype, float plane);

// ----------------------------------------------------------------------

struct mrc_io * ggcm_diag_lib_create_mrc_io(MPI_Comm comm, const char *run, const char *outputmode,
					    int outtype, float sheet, int rank_diagsrv);
void ggcm_diag_lib_write_openggcm_attrs(struct mrc_io *io, const char *time_str);
void ggcm_diag_lib_make_time_string(char s[80], float time, double dacttime);

// FIXME, this is specific to ggcm_mhd_diag "c" and shouldn't be here, but really,
// all of this hackiness should be resolved in some better way

// ======================================================================
// ggcm_mhd_diag_c

#define MAX_PLANES (10)

struct ggcm_mhd_diag_c {
  // parameters
  char *run;
  char *fields;
  char *outplanex;
  char *outplaney;
  char *outplanez;
  float planes[3][MAX_PLANES];
  int nr_planes[3];

  // state
  int rank_diagsrv;
  list_t mrc_io_list;

  // hacky indirection for OpenGGCM specific stuff
  void (*make_time_string)(char *time_str, float time, double dacttime);
  int  (*find_diagsrv)(struct ggcm_mhd_diag *diag);
  void (*run_hack)(struct ggcm_mhd_diag *diag, bool *output3d, bool *output2d,
		   int *itdia3d, int *itdia2d);
};

#define ggcm_mhd_diag_c(diag) mrc_to_subobj(diag, struct ggcm_mhd_diag_c)

#endif
