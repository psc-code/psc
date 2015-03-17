
#ifndef MRC_MAT_H
#define MRC_MAT_H

#include <mrc.h>
#include <mrc_domain.h>

// ======================================================================
// mrc_mat

MRC_CLASS_DECLARE(mrc_mat, struct mrc_mat);

void mrc_mat_assemble(struct mrc_mat *mat);
void mrc_mat_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val);
void mrc_mat_apply(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x);
void mrc_mat_apply_add(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x);
void mrc_mat_apply_in_place(struct mrc_mat *mat, struct mrc_vec *x);
void mrc_mat_apply_general(struct mrc_vec *z, double alpha,
			   struct mrc_mat *mat, struct mrc_vec *x,
			   double beta, struct mrc_vec *y);
void mrc_mat_print(struct mrc_mat *mat);

//======================================================================
// THE STORY OF THE FOLLOWING PART OF THIS FILE
//
// The code mrc-v3 had petsc integration before the rest of libmrc. It also had wonderful
// things like Multi-block domains and arbitrary metric coordinates. When it came time
// to bring it over, some of the things fit nicely, for example the multiblock domain just
// became a subclass of the domain. Other things, however, didn't. mrc-v3 needed matrix 
// support for it's use of the petsc implicit timestepper methods, so it has this simple framework 
// to wrap petsc matrices. I haven't gone beyond that, but someone should at some point.

#ifdef HAVE_PETSC

// ======================================================================
// matrix prealloc/create

#include <petscmat.h>

struct mat_create_ctx {
  int prealloc;
  int *dnz, *onz, rstart, rend, start, end, M, N;
};

int __MatCreate(MPI_Comm comm, int m, int n, int M, int N, Mat *mat,
		struct mat_create_ctx *mc);
int __MatSetValue(Mat M, int im, int in, PetscScalar v, int mode,
		  struct mat_create_ctx *mc);
int __MatInsertValue(Mat M, int im, int in, PetscScalar v,
		     struct mat_create_ctx *mc);



// To the poor soul who has to implement an mrc_matrix class: I have no idea why there is a 
// a public (mrc) and private (__mrc) version of this function, and if they are different.
// obviously the private version should never be used outside of the class, but for some reason
// they are both used. 

// Right now they live in mrc_ddc_mb.c. This is obviously the wrong place for them, but they
// require some ddc internals. One day maybe the hooks needed by these functions will
// exist in the ddc, but that's not something I have time for now.

// They used to be called MB_MatSetValue and __MB_MatSetValue, if you decide to crawl back through
// the code

int mrc_matrix_set_value(struct mrc_domain *mb, Mat mat, int bs, int im, int ig, 
		   int b, int jm, int jx, int jy, int jz,
		   double val);

int __mrc_matrix_set_value(struct mrc_domain *mb, Mat mat, int bs, int im, int ig, 
		     int b, int jm, int jx, int jy, int jz,
		     double val, struct mat_create_ctx *mc);

// This lives in ddc too, but it replaces a macro __I3g that was in the domain header
int mrc_matrix_find_global_index(struct mrc_domain *domain, 
				 int lpatch, 
				 int jx, int jy, int jz);

#endif

#endif
