
#include "mrc_mat_private.h"
#include "mrc_fld_as_double.h"
#include "mrc_vec.h"

#include "mrc_decomposition_private.h"

#include <stdlib.h>

// ======================================================================
// mrc_mat "csr_slow_mpi"

struct mrc_mat_csr_slow_mpi {
  struct mrc_mat *A;
  struct mrc_mat *B;

  struct mrc_decomposition *dc_row;
  struct mrc_decomposition *dc_col;

  int *rev_col_map; // FIXME, should go away, but still used for testing

  int n_recvs;
  int *recv_len;
  int *recv_src;
  struct mrc_vec *x_nl; // non-local compacted part of x that we need on the local proc

  int n_sends;
  int *send_len;
  int *send_dst;
  int *send_map;
  mrc_fld_data_t *send_buf;

  MPI_Request *req;
  bool is_assembled;
};

#define mrc_mat_csr_slow_mpi(mat) mrc_to_subobj(mat, struct mrc_mat_csr_slow_mpi)

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_create

static void
mrc_mat_csr_slow_mpi_create(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  sub->dc_row = mrc_decomposition_create(mrc_mat_comm(mat));
  sub->dc_col = mrc_decomposition_create(mrc_mat_comm(mat));

  sub->A = mrc_mat_create(MPI_COMM_SELF);
  mrc_mat_set_type(sub->A, "csr_slow");

  sub->B = mrc_mat_create(MPI_COMM_SELF);
  mrc_mat_set_type(sub->B, "csr_slow");
  
  sub->is_assembled = false;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_setup

static void
mrc_mat_csr_slow_mpi_setup(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  sub->dc_row->n = mat->m;
  sub->dc_col->n = mat->n;

  mrc_decomposition_setup(sub->dc_row);
  mrc_decomposition_setup(sub->dc_col);

  // A is the diagonal block, so use local sizes
  mrc_mat_set_param_int(sub->A, "m", sub->dc_row->n);
  mrc_mat_set_param_int(sub->A, "n", sub->dc_col->n);
  mrc_mat_setup(sub->A);

  // B is the off diagonal block, so # of rows = local # of rows,
  // but number of columns for now we set to the global number of columns
  // (even though local columns will never be inserted here, but
  // rather into A
  mrc_mat_set_param_int(sub->B, "m", sub->dc_row->n);
  mrc_mat_set_param_int(sub->B, "n", sub->dc_col->N);
  mrc_mat_setup(sub->B);

  mrc_mat_setup_super(mat);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_destroy

static void
mrc_mat_csr_slow_mpi_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  mrc_mat_destroy(sub->A);
  mrc_mat_destroy(sub->B);

  mrc_decomposition_destroy(sub->dc_row);
  mrc_decomposition_destroy(sub->dc_col);

  free(sub->recv_len);
  free(sub->recv_src);
  mrc_vec_destroy(sub->x_nl);

  free(sub->send_len);
  free(sub->send_dst);
  free(sub->send_map);
  free(sub->send_buf);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_add_value
//
// WARNING, all elements for any given row must be added contiguously!

static void
mrc_mat_csr_slow_mpi_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  assert(!sub->is_assembled);

  row_idx = mrc_decomposition_global_to_local(sub->dc_row, row_idx);

  if (mrc_decomposition_is_local(sub->dc_col, col_idx)) {
    col_idx = mrc_decomposition_global_to_local(sub->dc_col, col_idx);
    mrc_mat_add_value(sub->A, row_idx, col_idx, val);
  } else {
    mrc_mat_add_value(sub->B, row_idx, col_idx, val);
  }
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_assemble

static void
mrc_mat_csr_slow_mpi_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  assert(!sub->is_assembled);
  
  mrc_mat_assemble(sub->A);
  mrc_mat_assemble(sub->B);

  int rank, size;
  MPI_Comm_rank(mrc_mat_comm(mat), &rank);
  MPI_Comm_size(mrc_mat_comm(mat), &size);

  // create full ownership info on each proc
  int *col_offs_by_rank = calloc(size, sizeof(col_offs_by_rank));
  MPI_Allgather(&mat->n, 1, MPI_INT, col_offs_by_rank, 1, MPI_INT,
                mrc_mat_comm(mat));
  for (int r = 1; r < size; r++) {
    col_offs_by_rank[r] += col_offs_by_rank[r-1];
  }

  // set up map that maps non-zero column indices in B to a compacted index
  struct mrc_mat_csr_slow *sub_B = mrc_mat_csr_slow(sub->B);

  int N_cols = sub->dc_col->N;
  int *col_map = malloc(N_cols * sizeof(*col_map));
  for (int col = 0; col < N_cols; col++) {
    col_map[col] = -1;
  }

  // map each column with non-zero entries to an unique (per rank) entry in
  // the compact counting, starting at zero
  int *col_map_cnt_by_rank = calloc(size, sizeof(*col_map_cnt_by_rank));
  int col_map_cnt = 0;
  
  // FIXME: csr_slow specific
  for (int row=0; row < sub_B->nr_rows; row++) {
    for (int i=sub_B->rows[row]; i < sub_B->rows[row + 1]; i++) {
      int col = sub_B->cols[i];
      int r = mrc_decomposition_find_rank(sub->dc_col, col);
      if (col_map[col] == -1) {
        col_map[col] = col_map_cnt_by_rank[r]++;
        col_map_cnt++;
      }      
    }
  }

  // make col_map globally unique by adding the corresponding
  // per-rank offsets
  for (int r = 1; r < size; r++) {
    col_map_cnt_by_rank[r] += col_map_cnt_by_rank[r-1];
  }

  for (int col = 0; col < N_cols; col++) {
    if (col_map[col] >= 0) {
      int r = mrc_decomposition_find_rank(sub->dc_col, col);
      if (r > 0) {
        col_map[col] += col_map_cnt_by_rank[r-1];
      }
    }
  }

  free(col_map_cnt_by_rank);

  /* for (int col = 0; col < sub->N; col++) { */
  /*   if (col_map[col] >= 0) { */
  /*     mprintf("col_map %d -> %d\n", col, col_map[col]); */
  /*   } */
  /* } */

  // set up reverse map, mapping compacted indices back to original
  // global column indices
  sub->rev_col_map = malloc(col_map_cnt * sizeof(*sub->rev_col_map));

  for (int col = 0; col < N_cols; col++) {
    assert(col_map[col] < col_map_cnt);
    if (col_map[col] >= 0) {
      sub->rev_col_map[col_map[col]] = col;
    }
  }
  /* for (int i = 0; i < col_map_cnt; i++) { */
  /*   mprintf("rev map %d -> %d\n", i, sub->rev_col_map[i]); */
  /* } */

  // update column indices in B to refer to compacted indices
  // FIXME: csr_slow specific
  for (int row=0; row < sub_B->nr_rows; row++) {
    for (int i=sub_B->rows[row]; i < sub_B->rows[row + 1]; i++) {
      int col_idx = sub_B->cols[i];
      sub_B->cols[i] = col_map[col_idx];
    }
  }  
  sub->B->n = col_map_cnt;  // FIXME: this seems hacky, but doesn't break anything?
  free(col_map);

  // for each rank, find how many columns we need to receive from that rank
  sub->n_recvs = 0;
  int *recv_cnt_by_rank = calloc(size, sizeof(*recv_cnt_by_rank));
  for (int i = 0; i < col_map_cnt; i++) {
    int col = sub->rev_col_map[i];
    int r = mrc_decomposition_find_rank(sub->dc_col, col);
    if (recv_cnt_by_rank[r]++ == 0) {
      sub->n_recvs++;
    }
  }

  free(col_offs_by_rank);

  // we shouldn't have any columns to be received from our own rank,
  // because those are put into sub->A in the first place
  assert(recv_cnt_by_rank[rank] == 0);

  // compact recv_cnt_by_rank[]
  sub->recv_len = calloc(sub->n_recvs, sizeof(*sub->recv_len));
  sub->recv_src = calloc(sub->n_recvs, sizeof(*sub->recv_src));
  for (int n = 0, r = 0; r < size; r++) {
    if (recv_cnt_by_rank[r] > 0) {
      sub->recv_len[n] = recv_cnt_by_rank[r];
      sub->recv_src[n] = r;
      n++;
    }
  }

  free(recv_cnt_by_rank);

  // find out how many procs this proc needs to send messages to
  // (by finding this info for all procs first, then picking the
  // current rank's value)
  int *recv_flg_by_rank = calloc(size, sizeof(*recv_flg_by_rank));
  for (int n = 0; n < sub->n_recvs; n++) {
    recv_flg_by_rank[sub->recv_src[n]] = 1;
  }

  int *send_flg_by_rank = calloc(size, sizeof(*send_flg_by_rank));
  MPI_Allreduce(recv_flg_by_rank, send_flg_by_rank, size, MPI_INT, MPI_SUM,
                mrc_mat_comm(mat));
  sub->n_sends = send_flg_by_rank[rank];

  free(recv_flg_by_rank);
  free(send_flg_by_rank);

  // find compacted info on sends from the local proc,
  // by communication from those procs that expect to receive
  // those sends
  sub->send_len = calloc(sub->n_sends, sizeof(*sub->send_len));
  sub->send_dst = calloc(sub->n_sends, sizeof(*sub->send_dst));
  sub->req = calloc(sub->n_sends + sub->n_recvs, sizeof(*sub->req));
  mprintf("%p alloc sub->req %p %d %d\n", mat, sub->req, sub->n_sends, sub->n_recvs);
  for (int n = 0; n < sub->n_sends; n++) {
    MPI_Irecv(&sub->send_len[n], 1, MPI_INT, MPI_ANY_SOURCE, 0, mrc_mat_comm(mat),
              &sub->req[n]);
  }

  for (int n = 0; n < sub->n_recvs; n++) {
    MPI_Isend(&sub->recv_len[n], 1, MPI_INT, sub->recv_src[n], 0, mrc_mat_comm(mat),
              &sub->req[sub->n_sends + n]);
  }

  MPI_Status *status = calloc(sub->n_sends + sub->n_recvs, sizeof(*status));
  MHERE;
  MPI_Waitall(sub->n_sends + sub->n_recvs, sub->req, status);
  MHERE;
  for (int n = 0; n < sub->n_sends; n++) {
    sub->send_dst[n] = status[n].MPI_SOURCE;
  }
  free(status);

  /* for (int n = 0; n < sub->n_recvs; n++) { */
  /*   mprintf("recv_len[%d] = %d src = %d\n", n, sub->recv_len[n], sub->recv_src[n]); */
  /* } */
  /* for (int n = 0; n < sub->n_sends; n++) { */
  /*   mprintf("send_len[%d] = %d dst = %d\n", n, sub->send_len[n], sub->send_dst[n]); */
  /* } */

  // prepare actual send / recv buffers
  int send_buf_size = 0;
  for (int n = 0; n < sub->n_sends; n++) {
    send_buf_size += sub->send_len[n];
  }

  // now create a map on how to fill the send buffers
  sub->send_map = calloc(send_buf_size, sizeof(*sub->send_map));
  int *p = sub->send_map;
  for (int n = 0; n < sub->n_sends; n++) {
    MPI_Irecv(p, sub->send_len[n], MPI_INT, sub->send_dst[n], 1, mrc_mat_comm(mat),
              &sub->req[n]);
    p += sub->send_len[n];
  }
  p = sub->rev_col_map;
  for (int n = 0; n < sub->n_recvs; n++) {
    MPI_Isend(p, sub->recv_len[n], MPI_INT, sub->recv_src[n], 1, mrc_mat_comm(mat),
              &sub->req[sub->n_sends + n]);
    p += sub->recv_len[n];
  }

  MHERE;
  MPI_Waitall(sub->n_sends + sub->n_recvs, sub->req, MPI_STATUSES_IGNORE);
  MHERE;

  p = sub->send_map;
  for (int n = 0; n < sub->n_sends; n++) {
    for (int i = 0; i < sub->send_len[n]; i++) {
      p[i] = mrc_decomposition_global_to_local(sub->dc_col, p[i]);
      assert(p[i] >= 0 && p[i] < sub->dc_col->n);
      /* mprintf("map send %d: [%d] <- %d\n", n, i, p[i]); */
    }
    p += sub->send_len[n];
  }

  // compacted non-local part of right hand side vector x
  sub->x_nl = mrc_vec_create(mrc_mat_comm(mat));
  mrc_vec_set_type(sub->x_nl, FLD_TYPE);
  mrc_vec_set_param_int(sub->x_nl, "len", col_map_cnt);
  mrc_vec_setup(sub->x_nl);

  sub->send_buf = calloc(send_buf_size, sizeof(*sub->send_buf));
  
  sub->is_assembled = true;
}

static void
mrc_mat_csr_slow_mpi_gather_xc(struct mrc_mat *mat, struct mrc_vec *x, struct mrc_vec *x_nl)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  assert(sub->is_assembled);

  mrc_fld_data_t *x_arr = mrc_vec_get_array(x);
  mrc_fld_data_t *x_nl_arr = mrc_vec_get_array(x_nl);

#if 0
  // actual communication
  mrc_fld_data_t *pp = sub->send_buf;
  for (int n = 0; n < sub->n_sends; n++) {
    for (int i = 0; i < sub->send_len[n]; i++) {
      pp[i] = x_arr[sub->send_map[i]];
    }
    MPI_Isend(pp, sub->send_len[n], MPI_MRC_FLD_DATA_T, sub->send_dst[n], 2,
              mrc_mat_comm(mat), &sub->req[n]);
    pp += sub->send_len[n];
  }

  pp = x_nl_arr;
  for (int n = 0; n < sub->n_recvs; n++) {
    MPI_Irecv(pp, sub->recv_len[n], MPI_MRC_FLD_DATA_T, sub->recv_src[n], 2,
              mrc_mat_comm(mat), &sub->req[n]);
    pp += sub->recv_len[n];
  }

  MHERE;
  mprintf("%p sub->req %p %d %d\n", mat, sub->req, sub->n_sends, sub->n_recvs);
  MPI_Waitall(sub->n_sends + sub->n_recvs, sub->req, MPI_STATUSES_IGNORE);
  MHERE;
#endif

  // old
  // xg field -- FIXME, will be unneeded when we do proper communication
  struct mrc_vec *xg = mrc_vec_create(MPI_COMM_SELF);
  mrc_vec_set_type(xg, FLD_TYPE);
  mrc_vec_set_param_int(xg, "len", sub->dc_row->N);
  mrc_vec_setup(xg);

  int len = mrc_vec_len(x);
  mrc_fld_data_t *xg_arr = mrc_vec_get_array(xg);
  MPI_Allgather(x_arr, len, MPI_MRC_FLD_DATA_T,
                xg_arr, len, MPI_MRC_FLD_DATA_T, mrc_mat_comm(mat));

  for (int i = 0; i < mrc_vec_len(x_nl); i++) {
    x_nl_arr[i] = xg_arr[sub->rev_col_map[i]];
    //    assert(x_nl_arr[i] == xg_arr[sub->rev_col_map[i]]);
  }
  mrc_vec_put_array(xg, xg_arr);

  mrc_vec_destroy(xg);

  mrc_vec_put_array(x, x_arr);
  mrc_vec_put_array(x_nl, x_nl_arr);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_apply
// y = mat * x

static void
mrc_mat_csr_slow_mpi_apply(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  mrc_mat_apply(y, sub->A, x);
  mrc_mat_csr_slow_mpi_gather_xc(mat, x, sub->x_nl);
  mrc_mat_apply_add(y, sub->B, sub->x_nl);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_apply_in_place
// x = mat * x

static void
mrc_mat_csr_slow_mpi_apply_in_place(struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  mrc_mat_apply(x, sub->A, x);
  mrc_mat_csr_slow_mpi_gather_xc(mat, x, sub->x_nl);
  mrc_mat_apply_add(x, sub->B, sub->x_nl);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_apply_add
// y = mat * x + y

static void
mrc_mat_csr_slow_mpi_apply_add(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  mrc_mat_apply_add(y, sub->A, x);
  mrc_mat_csr_slow_mpi_gather_xc(mat, x, sub->x_nl);
  mrc_mat_apply_add(y, sub->B, sub->x_nl);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_apply_general
// z = alpha * mat * x + beta * y

static void
mrc_mat_csr_slow_mpi_apply_general(struct mrc_vec *z, double alpha,
                                   struct mrc_mat *mat, struct mrc_vec *x,
                                   double beta, struct mrc_vec *y)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  mrc_mat_apply_general(z, alpha, sub->A, x, beta, y);
  mrc_mat_csr_slow_mpi_gather_xc(mat, x, sub->x_nl);
  mrc_mat_apply_general(z, alpha, sub->B, sub->x_nl, 1.0, z);
}


// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi_print

static void
mrc_mat_csr_slow_mpi_print(struct mrc_mat *mat)
{
  struct mrc_mat_csr_slow_mpi *sub = mrc_mat_csr_slow_mpi(mat);

  mprintf("csr_slow_mpi sub-matrix A:\n");
  mrc_mat_print(sub->A);
  mprintf("\n");

  mprintf("csr_slow_mpi sub-matrix B:\n");
  mrc_mat_print(sub->B);
  mprintf("\n");
}

// ----------------------------------------------------------------------
// mrc_mat_csr_slow_mpi description

#define VAR(x) (void *)offsetof(struct mrc_mat_csr_slow_mpi, x)
static struct param mrc_mat_csr_slow_mpi_descr[] = {
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_mat subclass "csr_slow_mpi"

struct mrc_mat_ops mrc_mat_csr_slow_mpi_ops = {
  .name                  = "csr_slow_mpi",
  .size                  = sizeof(struct mrc_mat_csr_slow_mpi),
  .param_descr           = mrc_mat_csr_slow_mpi_descr,
  .create                = mrc_mat_csr_slow_mpi_create,
  .setup                 = mrc_mat_csr_slow_mpi_setup,
  .destroy               = mrc_mat_csr_slow_mpi_destroy,
  .add_value             = mrc_mat_csr_slow_mpi_add_value,
  .assemble              = mrc_mat_csr_slow_mpi_assemble,
  .apply                 = mrc_mat_csr_slow_mpi_apply,
  .apply_in_place        = mrc_mat_csr_slow_mpi_apply_in_place,
  .apply_add             = mrc_mat_csr_slow_mpi_apply_add,
  .apply_general         = mrc_mat_csr_slow_mpi_apply_general,
  .print                 = mrc_mat_csr_slow_mpi_print,
};
