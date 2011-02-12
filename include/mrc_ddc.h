
#ifndef DDC_H
#define DDC_H

#include <mpi.h>

enum {
  DDC_BC_NONE,
  DDC_BC_PERIODIC,
};

struct mrc_ddc_params {
  MPI_Comm comm;
  MPI_Datatype mpi_type;
  int size_of_type;
  int max_n_fields;
  int n_proc[3]; // # procs in 3D grid
  int ilo[3], ihi[3]; // local domain (no ghosts)
  int ibn[3]; // # ghost points
  int bc[3]; // boundary condition
  void (*copy_to_buf)(int mb, int me, int ilo[3], int ihi[3], void *buf, void *ctx);
  void (*copy_from_buf)(int mb, int me, int ilo[3], int ihi[3], void *buf, void *ctx);
  void (*add_from_buf)(int mb, int me, int ilo[3], int ihi[3], void *buf, void *ctx);
};

struct mrc_ddc_sendrecv {
  int ilo[3], ihi[3];
  int rank_nei;
  int len;
  void *buf;
};

struct mrc_ddc_pattern {
  struct mrc_ddc_sendrecv send[27];
  struct mrc_ddc_sendrecv recv[27];
};

struct mrc_ddc {
  struct mrc_ddc_params prm;
  int rank, size;
  int proc[3]; // this proc's position in the 3D proc grid
  struct mrc_ddc_pattern add_ghosts;
  struct mrc_ddc_pattern fill_ghosts;
  MPI_Request send_reqs[27];
  MPI_Request recv_reqs[27];
};

struct mrc_ddc *mrc_ddc_create(struct mrc_ddc_params *prm);
void mrc_ddc_add_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx);
void mrc_ddc_fill_ghosts(struct mrc_ddc *ddc, int mb, int me, void *ctx);

int mrc_ddc_get_rank_nei(struct mrc_ddc *ddc, int dir[3]);

#define MRC_DDC_BUF3(buf,m, ix,iy,iz)		\
  (buf[(((m) * (ihi[2] - ilo[2]) +		\
	 iz - ilo[2]) * (ihi[1] - ilo[1]) +	\
	iy - ilo[1]) * (ihi[0] - ilo[0]) +	\
       ix - ilo[0]])

#endif

static inline int
mrc_ddc_dir2idx(int dir[3])
{
  return ((dir[2] + 1) * 3 + dir[1] + 1) * 3 + dir[0] + 1;
}

