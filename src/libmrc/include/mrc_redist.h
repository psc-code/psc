
#ifndef MRC_REDIST_H
#define MRC_REDIST_H

#include <mrc.h>
#include <mrc_fld.h>

BEGIN_C_DECLS

struct mrc_redist_block {
  int ilo[3]; // intersection low
  int ihi[3]; // intersection high
  int p; // local patch
};

struct mrc_redist_peer {
  int rank;
  struct mrc_redist_block *blocks_begin, *blocks_end;
  // off and buf_size are redundant with disps
  int off;
  int buf_size;
};

struct mrc_redist_write_recv {
  struct mrc_redist_peer *peers_begin, *peers_end;
  MPI_Request *reqs;
  void *buf;
  size_t buf_size;
  int *cnts;
  int *disps;

  int n_recv_patches;
  struct mrc_redist_block *recv_patches;
};

struct mrc_redist_write_send {
  struct mrc_redist_peer *peers_begin, *peers_end;
  MPI_Request *reqs;
  void *buf;
  size_t buf_size;
  int *cnts;
  int *disps;
};

struct mrc_redist {
  struct mrc_domain *domain;
  MPI_Comm comm;
  int rank;
  int size;
  
  MPI_Comm comm_writers;
  int *writer_ranks;
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

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

static inline bool
find_intersection(int *ilo, int *ihi, const int *ib1, const int *im1,
		  const int *ib2, const int *im2)
{
  bool has_intersection = true;
  for (int d = 0; d < 3; d++) {
    ilo[d] = MAX(ib1[d], ib2[d]);
    ihi[d] = MIN(ib1[d] + im1[d], ib2[d] + im2[d]);
    if (ihi[d] - ilo[d] <= 0) {
      has_intersection = false;
    }
  }
  return has_intersection;
}

END_C_DECLS

#endif
