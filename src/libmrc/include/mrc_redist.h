
#ifndef MRC_REDIST_H
#define MRC_REDIST_H

BEGIN_C_DECLS

struct mrc_redist_block {
  int ilo[3]; // intersection low
  int ihi[3]; // intersection high
};

struct mrc_redist_peer {
  int rank;
  struct mrc_redist_block *begin;
  struct mrc_redist_block *end;
  void *buf;
  size_t buf_size;
};

struct mrc_redist_write_recv {
  int n_peers;
  struct mrc_redist_peer *peers;
  MPI_Request *reqs;

  int n_recv_patches;
  struct mrc_redist_block *recv_patches;
};

struct mrc_redist_write_send {
  void **bufs; // one for each writer
  int *buf_sizes;
  MPI_Request *reqs;
};

struct mrc_redist {
  struct mrc_domain *domain;
  MPI_Comm comm;
  int rank;
  int size;
  
  MPI_Comm comm_writers;
  int *writers;
  int nr_writers;
  int is_writer;

  int slab_dims[3];
  int slab_offs[3];

  int slow_dim;
  int slow_indices_per_writer;
  int slow_indices_rmndr;

  struct mrc_redist_write_send write_send;
  struct mrc_redist_write_recv write_recv;
};

void mrc_redist_init(struct mrc_redist *redist, struct mrc_domain *domain,
		     int slab_offs[3], int slab_dims[3], int nr_writers);
void mrc_redist_destroy(struct mrc_redist *redist);
struct mrc_ndarray *mrc_redist_get_ndarray(struct mrc_redist *redist, struct mrc_fld *m3);
void mrc_redist_put_ndarray(struct mrc_redist *redist, struct mrc_ndarray *nd);
void mrc_redist_run(struct mrc_redist *redist, struct mrc_ndarray *nd,
		    struct mrc_fld *m3, int m);

END_C_DECLS

#endif
