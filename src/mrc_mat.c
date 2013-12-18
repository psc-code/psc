
#include <mrc_mat.h>
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
