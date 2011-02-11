
#ifndef DDC_H
#define DDC_H

#include <mrc_common.h>

#include <mpi.h>

struct mrc_ddc;

struct mrc_ddc_params {
  int n[3];       // size of local domain (w/o ghosts)
  int np[3];      // number of processors in each crd direction
  int bc[3];      // boundary condition in each crd direction (BC_NONE,BC_PERIODIC)
  int bnd;        // number of ghost points
  int max_n_comp; // maximum # of components in fields
};

struct mrc_ddc_ops {
  void (*copy_to_buf)(struct mrc_ddc_params *ddc_params, float *x, float *buf,
		      const int ib[3], const int ie[3], int mb, int me);
  void (*copy_from_buf)(struct mrc_ddc_params *ddc_params, float *x, float *buf,
			const int ib[3], const int ie[3], int mb, int me);
};

struct mrc_ddc *mrc_ddc_create(MPI_Comm, struct mrc_ddc_params *params,
			       struct mrc_ddc_ops *ops);
void mrc_ddc_destroy(struct mrc_ddc *ddc);
void mrc_ddc_fill_ghosts_begin(struct mrc_ddc *ddc, float *x, int mb, int me);
void mrc_ddc_fill_ghosts_end(struct mrc_ddc *ddc, float *x, int mb, int me);
void mrc_ddc_fill_ghosts(struct mrc_ddc *ddc, float *x, int mb, int me);

// for implementing copy_{to,from}_buf()

#define MRC_DDC_BUF3(buf,m, ix,iy,iz)					\
  ((buf)[((((iz)-ib[2]) * (ie[1]-ib[1]) + ((iy)-ib[1])) * (ie[0]-ib[0]) + ((ix)-ib[0])) * (me-mb) + ((m)-mb)])

// for standard Fortran layout fields (component: slow idx)
extern struct mrc_ddc_ops mrc_ddc_ops_fortran;

#endif
