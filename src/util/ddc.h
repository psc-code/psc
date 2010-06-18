
#ifndef DDC_H
#define DDC_H

#include <mpi.h>

typedef double ddc_real; // FIXME

enum {
  DDC_BC_NONE,
  DDC_BC_PERIODIC,
};

struct ddc_params {
  MPI_Comm comm;
  int n_proc[3]; // # procs in 3D grid
  int ilo[3], ihi[3]; // local domain (no ghosts)
  int ibn[3]; // # ghost points
  int bc[3]; // boundary condition
  void (*copy_to_buf)(int m, int ilo[3], int ihi[3], ddc_real *buf);
  void (*add_from_buf)(int m, int ilo[3], int ihi[3], ddc_real *buf);
};

struct ddc_send {
  int ilo[3], ihi[3];
  int rank_nei;
  int len;
  ddc_real *buf;
};

struct ddc_recv {
  int ilo[3], ihi[3];
  int rank_nei;
  int len;
  ddc_real *buf;
};

struct ddc_subdomain {
  struct ddc_params prm;
  int rank, size;
  int proc[3]; // this proc's position in the 3D proc grid
  struct ddc_send send[27];
  struct ddc_recv recv[27];
  MPI_Request send_reqs[27];
  MPI_Request recv_reqs[27];
};

struct ddc_subdomain *ddc_create(struct ddc_params *prm);
void ddc_add_ghosts(struct ddc_subdomain *ddc, int m);

#define DDC_BUF(buf, ix,iy,iz)			\
  (buf[((iz - ilo[2]) * (ihi[1] - ilo[1]) +	\
	iy - ilo[1]) * (ihi[0] - ilo[0]) +	\
       ix - ilo[0]])

#endif

