
#ifndef MRC_DDC_PRIVATE_H
#define MRC_DDC_PRIVATE_H

#include "mrc_ddc.h"

#define N_DIR (27)

struct mrc_ddc_sendrecv {
  int ilo[3], ihi[3];
  int nei_rank;
  int nei_patch;
  int len;
  void *buf;
};

struct mrc_ddc {
  struct mrc_obj obj;
  // parameters
  int size_of_type;
  int max_n_fields;
  int ibn[3]; // # ghost points

  int rank, size;
  struct mrc_ddc_funcs *funcs;
  MPI_Datatype mpi_type;
};

struct mrc_ddc_ops {
  MRC_SUBCLASS_OPS(struct mrc_ddc);
  void (*set_domain)(struct mrc_ddc *ddc, struct mrc_domain *domain);
  struct mrc_domain *(*get_domain)(struct mrc_ddc *ddc);
  void (*fill_ghosts_fld)(struct mrc_ddc *ddc, int mb, int me,
			  struct mrc_fld *fld);
  // FIXME: Needed for MB and nothing else!
  void (*global_to_local_fld)(struct mrc_ddc *ddc, struct mrc_fld *gfld,
			      struct mrc_fld *lfld);
  void (*fill_ghost_edges_fld)(struct mrc_ddc *ddc, int mb, int me,
			  struct mrc_fld *fld);  
  // OBSOLETE
  void (*fill_ghosts)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*fill_ghosts_begin)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*fill_ghosts_end)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*fill_ghosts_local)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
  void (*add_ghosts)(struct mrc_ddc *ddc, int mb, int me, void *ctx);
};

extern struct mrc_ddc_ops mrc_ddc_simple_ops;
extern struct mrc_ddc_ops mrc_ddc_multi_ops;
extern struct mrc_ddc_ops mrc_ddc_amr_ops;
extern struct mrc_ddc_ops mrc_ddc_mb_ops;
#endif
