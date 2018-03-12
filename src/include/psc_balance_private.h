
#ifndef PSC_BALANCE_PRIVATE_H
#define PSC_BALANCE_PRIVATE_H

#include <psc_balance.h>

struct psc_balance {
  struct mrc_obj obj;
  int every;
  double factor_fields;
  bool print_loads;
  bool write_loads;
};

struct communicate_ctx;

struct psc_balance_ops {
  MRC_SUBCLASS_OPS(struct psc_balance);
  const char *mprts_type;
  const char *mflds_type;
};

#define psc_balance_ops(bal) ((struct psc_balance_ops *)(bal->obj.ops))

struct send_info {
  int rank;
  int patch;
};

struct recv_info {
  int rank;
  int patch;
};

struct by_ri {
  int rank;
  int nr_patches;
  int *pi_to_patch;
};

struct communicate_ctx {
  MPI_Comm comm;
  int mpi_rank;
  int mpi_size;
  int nr_patches_old;
  int nr_patches_new;
  struct send_info *send_info; // by old patch on this proc
  struct recv_info *recv_info; // by new patch on this proc

  int *send_rank_to_ri; // map from send target rank to contiguous "rank index"
  int *recv_rank_to_ri; // map from recv source rank to contiguous "rank index"

  int nr_send_ranks;
  struct by_ri *send_by_ri;

  int nr_recv_ranks;
  struct by_ri *recv_by_ri;
};

#endif
