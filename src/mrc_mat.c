
#include <mrc_mat.h>

#include <stdlib.h>

// ======================================================================
// mrc_mat "mcsr"

// ----------------------------------------------------------------------
// mrc_mat_setup

void
_mrc_mat_setup(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = &mat->sub;

  sub->nr_rows_alloced = 1000;
  sub->nr_entries_alloced = 2000;

  sub->rows = calloc(sub->nr_rows_alloced, sizeof(*sub->rows));
  sub->entries = calloc(sub->nr_entries_alloced, sizeof(*sub->entries));

  sub->nr_entries = 0;
  sub->nr_rows = 0;
}

// ----------------------------------------------------------------------
// mrc_mat_destroy

void
_mrc_mat_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = &mat->sub;

  free(sub->rows);
  free(sub->entries);
}

// ----------------------------------------------------------------------
// mrc_mat_add_value

void
mrc_mat_add_value(struct mrc_mat *mat, int row_idx, int col_idx, float val)
{
  struct mrc_mat_mcsr *sub = &mat->sub;

  if (sub->nr_rows == 0 ||
      sub->rows[sub->nr_rows - 1].idx != row_idx) {
    // start new row
    if (sub->nr_rows >= sub->nr_rows_alloced - 1) {
      sub->nr_rows_alloced *= 2;
      sub->rows = realloc(sub->rows, sub->nr_rows_alloced * sizeof(*sub->rows));
    }
    sub->rows[sub->nr_rows].idx = row_idx;
    sub->rows[sub->nr_rows].first_entry = sub->nr_entries;
    sub->nr_rows++;
  }

  // if we already have an entry for this column in the current row, just add to it
  for (int i = sub->rows[sub->nr_rows - 1].first_entry; i < sub->nr_entries; i++) {
    if (sub->entries[i].idx == col_idx) {
      sub->entries[i].val += val;
      return;
    }
  }

  // otherwise, need to append a new entry
  if (sub->nr_entries >= sub->nr_entries_alloced) {
    sub->nr_entries_alloced *= 2;
    sub->entries = realloc(sub->entries, sub->nr_entries_alloced * sizeof(*sub->entries));
  }
  sub->entries[sub->nr_entries].idx = col_idx;
  sub->entries[sub->nr_entries].val = val;
  sub->nr_entries++;
}

// ----------------------------------------------------------------------
// mrc_mat_assemble

void
mrc_mat_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr *sub = &mat->sub;

  sub->rows[sub->nr_rows].first_entry = sub->nr_entries;
  mprintf("nr_rows %d nr_entries %d\n", sub->nr_rows, sub->nr_entries);
}

// ----------------------------------------------------------------------
// mrc_mat_apply

void
mrc_mat_apply(struct mrc_mat *mat, struct mrc_fld *fld)
{
  struct mrc_mat_mcsr *sub = &mat->sub;

  float *arr = fld->_arr;
    
  for (int row = 0; row < sub->nr_rows; row++) {
    int row_idx = sub->rows[row].idx;
    float sum = 0.;
    for (int entry = sub->rows[row].first_entry;
	 entry < sub->rows[row + 1].first_entry; entry++) {
      int col_idx = sub->entries[entry].idx;
      float val = sub->entries[entry].val;
      sum += val * arr[col_idx];
    }
    arr[row_idx] = sum;
  }
}

// ----------------------------------------------------------------------
// mrc_mat class description

struct mrc_class_mrc_mat mrc_class_mrc_mat = {
  .name         = "mrc_mat",
  .size         = sizeof(struct mrc_mat),
  .setup        = _mrc_mat_setup,
  .destroy      = _mrc_mat_destroy,
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
