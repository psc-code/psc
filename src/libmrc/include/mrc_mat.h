
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


#endif

#endif
