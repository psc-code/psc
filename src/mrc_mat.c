
#include <mrc_mat_private.h>

#include <stdlib.h>

#define mrc_mat_ops(mat) ((struct mrc_mat_ops *) mat->obj.ops)

// ----------------------------------------------------------------------
// mrc_mat_add_value

void
mrc_mat_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->add_value);
  ops->add_value(mat, row_idx, col_idx, val);
}

// ----------------------------------------------------------------------
// mrc_mat_assemble

void
mrc_mat_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->assemble);
  ops->assemble(mat);
}

// ----------------------------------------------------------------------
// mrc_mat_apply

void
mrc_mat_apply(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->apply);
  ops->apply(y, mat, x);
}

// ----------------------------------------------------------------------
// mrc_mat_apply_add

void
mrc_mat_apply_add(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->apply_add);
  ops->apply_add(y, mat, x);
}

// ----------------------------------------------------------------------
// mrc_mat_apply_in_place

void
mrc_mat_apply_in_place(struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->apply_in_place);
  ops->apply_in_place(mat, x);
}

// ----------------------------------------------------------------------
// mrc_mat_apply_general

void mrc_mat_apply_general(struct mrc_vec *z, double alpha,
			   struct mrc_mat *mat, struct mrc_vec *x,
			   double beta, struct mrc_vec *y)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->apply_general);
  ops->apply_general(z, alpha, mat, x, beta, y);
}

// ----------------------------------------------------------------------
// mrc_mat_print

void
mrc_mat_print(struct mrc_mat *mat)
{
  struct mrc_mat_ops *ops = mrc_mat_ops(mat);
  assert(ops->print);
  ops->print(mat);
}

// ----------------------------------------------------------------------
// mrc_mat_init

static void
mrc_mat_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_mat, &mrc_mat_csr_ops);
  mrc_class_register_subclass(&mrc_class_mrc_mat, &mrc_mat_csr_mpi_ops);
  mrc_class_register_subclass(&mrc_class_mrc_mat, &mrc_mat_mcsr_ops);
  mrc_class_register_subclass(&mrc_class_mrc_mat, &mrc_mat_mcsr_mpi_ops);
#ifdef HAVE_PETSC
  mrc_class_register_subclass(&mrc_class_mrc_mat, &mrc_mat_petsc_ops);
#endif
}

// ----------------------------------------------------------------------
// mrc_mat description

#define VAR(x) (void *)offsetof(struct mrc_mat, x)
static struct param mrc_mat_descr[] = {
  { "m"                 , VAR(m)                 , PARAM_INT(0) },
  { "n"                 , VAR(n)                 , PARAM_INT(0) },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_mat class

struct mrc_class_mrc_mat mrc_class_mrc_mat = {
  .name         = "mrc_mat",
  .size         = sizeof(struct mrc_mat),
  .param_descr  = mrc_mat_descr,
  .init         = mrc_mat_init,
};


// ======================================================================
// petsc-specific function that should be revisited eventually FIXME

#ifdef HAVE_PETSC

#define CE assert(ierr == 0)

int
__MatCreate(MPI_Comm comm, int m, int n, int M, int N, Mat *mat,
	    struct mat_create_ctx *ctx)
{
  int ierr;

  PetscFunctionBegin;
  if (ctx->prealloc == 0) {
    ierr = MatCreate(comm, mat); CE;
    ierr = MatSetSizes(*mat, m, n, M, N); CE;
    ierr = MatSetFromOptions(*mat); CE;

    ierr = MatGetSize(*mat, &M, &N); CE;
    ierr = MatGetLocalSize(*mat, &m, &n); CE;
    ctx->M = M; ctx->N = N;

    ierr = PetscMalloc(2*m*sizeof(int), &ctx->dnz); CE; 
    ctx->onz = ctx->dnz + m;
    ierr = PetscMemzero(ctx->dnz, 2*m*sizeof(int)); CE;
    ierr = MPI_Scan(&n, &ctx->end, 1, MPI_INT, MPI_SUM, comm); CE;
    ctx->start = ctx->end - n;
    ierr = MPI_Scan(&m, &ctx->rend, 1, MPI_INT, MPI_SUM, comm); CE;
    ctx->rstart = ctx->rend - m;
  } else {
    ierr = MatSeqAIJSetPreallocation(*mat, 0, ctx->dnz); CE;
    ierr = MatMPIAIJSetPreallocation(*mat, 0, ctx->dnz, 0, ctx->onz); CE;

    ierr = PetscFree(ctx->dnz); CE;
  }
  PetscFunctionReturn(0);
}

int
__MatSetValue(Mat M, int im, int in, PetscScalar v, int mode,
	      struct mat_create_ctx *ctx)
{
  int ierr, rank;
  MPI_Comm comm;

  PetscFunctionBegin;

  ierr = PetscObjectGetComm((PetscObject) M, &comm); CE;
  ierr = MPI_Comm_rank(comm, &rank); CE;

  if (ctx->prealloc == 0) {
    if (im < ctx->rstart || im >= ctx->rend) {
      SETERRQ3(MPI_COMM_WORLD, 1, "im %d out of bounds [%d:%d]", im, ctx->rstart, ctx->rend);
    }
    if (in < 0 || in >= ctx->N) {
      SETERRQ3(MPI_COMM_WORLD, 1, "in %d out of bounds [%d;%d]", in, 0, ctx->N);
    }
   
    if (in < ctx->start || in >= ctx->end) {
      ctx->onz[im - ctx->rstart]++;
    } else {
      ctx->dnz[im - ctx->rstart]++;
    }
  } else {
    ierr = MatSetValue(M, im, in, v, mode); CE;
  }

  PetscFunctionReturn(0);
}

int
__MatInsertValue(Mat M, int im, int in, PetscScalar v,
		 struct mat_create_ctx *ctx)
{
  return __MatSetValue(M, im, in, v, INSERT_VALUES, ctx);  
}

#endif
