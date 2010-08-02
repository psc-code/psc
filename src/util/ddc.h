
#ifndef DDC_H
#define DDC_H

#include <mpi.h>

enum {
  DDC_BC_NONE,
  DDC_BC_PERIODIC,
};

struct ddc_params {
  MPI_Comm comm;
  MPI_Datatype mpi_type;
  int size_of_type;
  int max_n_fields;
  int n_proc[3]; // # procs in 3D grid
  int ilo[3], ihi[3]; // local domain (no ghosts)
  int ibn[3]; // # ghost points
  int bc[3]; // boundary condition
  void (*copy_to_buf)(int m, int ilo[3], int ihi[3], void *buf);
  void (*copy_from_buf)(int m, int ilo[3], int ihi[3], void *buf);
  void (*add_from_buf)(int m, int ilo[3], int ihi[3], void *buf);
};

struct ddc_sendrecv {
  int ilo[3], ihi[3];
  int rank_nei;
  int len;
  void *buf;
};

struct ddc_pattern {
  struct ddc_sendrecv send[27];
  struct ddc_sendrecv recv[27];
};

struct ddc_subdomain {
  struct ddc_params prm;
  int rank, size;
  int proc[3]; // this proc's position in the 3D proc grid
  struct ddc_pattern add_ghosts;
  struct ddc_pattern fill_ghosts;
  MPI_Request send_reqs[27];
  MPI_Request recv_reqs[27];
};

struct ddc_subdomain *ddc_create(struct ddc_params *prm);
void ddc_add_ghosts(struct ddc_subdomain *ddc, int m);
void ddc_fill_ghosts(struct ddc_subdomain *ddc, int m);

int ddc_get_rank_nei(struct ddc_subdomain *ddc, int dir[3]);

#define DDC_BUF(buf, ix,iy,iz)			\
  (buf[((iz - ilo[2]) * (ihi[1] - ilo[1]) +	\
	iy - ilo[1]) * (ihi[0] - ilo[0]) +	\
       ix - ilo[0]])

#endif

static inline int
dir2idx(int dir[3])
{
  return ((dir[2] + 1) * 3 + dir[1] + 1) * 3 + dir[0] + 1;
}

