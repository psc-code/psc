
#include "mrc_mat_private.h"
#include "mrc_fld_as_double.h"
#include "mrc_ddc_private.h" 
#include <stdlib.h>

// ======================================================================
// mrc_mat "mcsr_mpi"

struct mrc_mat_mcsr_mpi {
  struct mrc_mat *A;
  struct mrc_mat *B;

  int M; // global number of rows
  int N; // global number of columns
  int m_off; // gobal index of first (0th) row that lives on this proc
  int n_off; // gobal index of first (0th) column that lives on this proc
};

#define mrc_mat_mcsr_mpi(mat) mrc_to_subobj(mat, struct mrc_mat_mcsr_mpi)

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_create

static void
mrc_mat_mcsr_mpi_create(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  sub->A = mrc_mat_create(MPI_COMM_SELF);
  mrc_mat_set_type(sub->A, "mcsr");

  sub->B = mrc_mat_create(MPI_COMM_SELF);
  mrc_mat_set_type(sub->B, "mcsr");
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_setup

static void
mrc_mat_mcsr_mpi_setup(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  MPI_Allreduce(&mat->m, &sub->M, 1, MPI_INT, MPI_SUM, mrc_mat_comm(mat));
  MPI_Allreduce(&mat->n, &sub->N, 1, MPI_INT, MPI_SUM, mrc_mat_comm(mat));
  MPI_Exscan(&mat->m, &sub->m_off, 1, MPI_INT, MPI_SUM, mrc_mat_comm(mat));
  MPI_Exscan(&mat->n, &sub->n_off, 1, MPI_INT, MPI_SUM, mrc_mat_comm(mat));

  mprintf("M = %d N = %d\n", sub->M, sub->N);
  mprintf("m_off = %d n_off = %d\n", sub->m_off, sub->n_off);

  // A is the diagonal block, so use local sizes
  mrc_mat_set_param_int(sub->A, "m", mat->m);
  mrc_mat_set_param_int(sub->A, "n", mat->n);
  mrc_mat_setup(sub->A);

  // B is the off diagonal block, so # of rows = local # of rows,
  // but number of columns for now we set to the global number of columns
  // (even though local columns will never be inserted here, but
  // rather into A
  mrc_mat_set_param_int(sub->B, "m", mat->m);
  mrc_mat_set_param_int(sub->B, "n", sub->N);
  mrc_mat_setup(sub->B);

  mrc_mat_setup_super(mat);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_destroy

static void
mrc_mat_mcsr_mpi_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  mrc_mat_destroy(sub->A);
  mrc_mat_destroy(sub->B);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_add_value
//
// WARNING, all elements for any given row must be added contiguously!

static void
mrc_mat_mcsr_mpi_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  row_idx -= sub->m_off;
  assert(row_idx >= 0 && row_idx < mat->m);
  
  if (col_idx >= sub->n_off && col_idx < sub->n_off + mat->n) {
    col_idx -= sub->n_off;
    mrc_mat_add_value(sub->A, row_idx, col_idx, val);
  } else {
    mrc_mat_add_value(sub->B, row_idx, col_idx, val);
  }
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_assemble

static void
mrc_mat_mcsr_mpi_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  mrc_mat_assemble(sub->A);
  mrc_mat_assemble(sub->B);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_apply

static void
mrc_mat_mcsr_mpi_apply(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);
  struct mrc_mat_mcsr *sub_B = mrc_mat_mcsr(sub->B);

  int *col_map = malloc(sub->N * sizeof(*col_map));
  for (int col = 0; col < sub->N; col++) {
    col_map[col] = -1;
  }
  int col_map_cnt = 0;
  for (int row = 0; row < sub_B->nr_rows; row++) {
    for (int entry = sub_B->rows[row].first_entry;
	 entry < sub_B->rows[row + 1].first_entry; entry++) {
      int col_idx = sub_B->entries[entry].idx;
      if (col_map[col_idx] == -1) {
	col_map[col_idx] = col_map_cnt++;
      }
    }
  }
  for (int col = 0; col < sub->N; col++) {
    mprintf("map %d -> %d\n", col, col_map[col]);
  }
  int *rev_col_map = malloc(col_map_cnt * sizeof(*rev_col_map));
  for (int col = 0; col < sub->N; col++) {
    rev_col_map[col_map[col]] = col;
  }
  for (int i = 0; i < col_map_cnt; i++) {
    mprintf("rev map %d -> %d\n", i, rev_col_map[i]);
  }

  struct mrc_fld *xg = mrc_fld_create(mrc_mat_comm(mat));
  mrc_fld_set_type(xg, FLD_TYPE);
  mrc_fld_set_param_int_array(xg, "dims", 1, (int[1]) { sub->M });
  mrc_fld_setup(xg);

  MPI_Allgather(x->_arr, x->_len, MPI_MRC_FLD_DATA_T,
		xg->_arr, x->_len, MPI_MRC_FLD_DATA_T, mrc_mat_comm(mat));

  struct mrc_fld *xc = mrc_fld_create(mrc_mat_comm(mat));
  mrc_fld_set_type(xc, FLD_TYPE);
  mrc_fld_set_param_int_array(xc, "dims", 1, (int[1]) { col_map_cnt });
  mrc_fld_setup(xc);

  for (int i = 0; i < col_map_cnt; i++) {
    MRC_D1(xc, i) = MRC_D1(xg, rev_col_map[i]);
    mprintf("xc[%d] = %g\n", i, MRC_D1(xc, i));
  }

  for (int row = 0; row < sub_B->nr_rows; row++) {
    for (int entry = sub_B->rows[row].first_entry;
	 entry < sub_B->rows[row + 1].first_entry; entry++) {
      int col_idx = sub_B->entries[entry].idx;
      sub_B->entries[entry].idx = col_map[col_idx];
    }
  }

  mrc_mat_apply(y, sub->A, x);
  mrc_mat_apply_add(y, sub->B, xc);

  mrc_fld_destroy(xg);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_print

static void
mrc_mat_mcsr_mpi_print(struct mrc_mat *mat)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  mprintf("mcsr_mpi sub-matrix A:\n");
  mrc_mat_print(sub->A);
  mprintf("\n");

  mprintf("mcsr_mpi sub-matrix B:\n");
  mrc_mat_print(sub->B);
  mprintf("\n");
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi description

#define VAR(x) (void *)offsetof(struct mrc_mat_mcsr_mpi, x)
static struct param mrc_mat_mcsr_mpi_descr[] = {
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_mat subclass "mcsr_mpi"

struct mrc_mat_ops mrc_mat_mcsr_mpi_ops = {
  .name                  = "mcsr_mpi",		
  .size                  = sizeof(struct mrc_mat_mcsr_mpi),
  .param_descr           = mrc_mat_mcsr_mpi_descr,
  .create                = mrc_mat_mcsr_mpi_create,
  .setup                 = mrc_mat_mcsr_mpi_setup,
  .destroy               = mrc_mat_mcsr_mpi_destroy,
  .add_value             = mrc_mat_mcsr_mpi_add_value,
  .assemble              = mrc_mat_mcsr_mpi_assemble,
  .apply                 = mrc_mat_mcsr_mpi_apply,
  .print                 = mrc_mat_mcsr_mpi_print,
  /* .apply_in_place        = mrc_mat_mcsr_mpi_apply_in_place, */
};

