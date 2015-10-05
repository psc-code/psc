
#include "mrc_mat_private.h"
#include "mrc_fld_as_double.h"
#include "mrc_vec.h"
#include "mrc_profile.h"

#include "mrc_decomposition_private.h"

#include <stdlib.h>

// ======================================================================
// mrc_mat "csr_mpi"

struct mrc_mat_csr_mpi {
  struct mrc_mat *A;
  struct mrc_mat *B;

  struct mrc_decomposition *dc_row;
  struct mrc_decomposition *dc_col;

  int *rev_col_map; // FIXME, should go away, but still used for testing

  int n_recvs;
  int *recv_len;
  int *recv_src;
  
  // non-local compacted part of x that we need on the local proc
  // during communication, this is set to NULL to indicate that it's
  // off limits
  struct mrc_vec *x_nl;

  int n_sends;
  int *send_len;
  int *send_dst;
  int *send_map;
  
  mrc_fld_data_t *send_buf;
  mrc_fld_data_t *_recv_buf;
  struct mrc_vec *_recv_vec;

  MPI_Request *req;
  bool is_assembled;
  
  bool verbose;
  bool do_profiling;
  int nr_initial_cols;  // how much initial space to allocate for each new row
  float growth_factor;  // fraction of current # of cols to add when a row needs to grow
};

#define mrc_mat_csr_mpi(mat) mrc_to_subobj(mat, struct mrc_mat_csr_mpi)

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_create

static void
mrc_mat_csr_mpi_create(struct mrc_mat *mat)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  sub->dc_row = mrc_decomposition_create(mrc_mat_comm(mat));
  sub->dc_col = mrc_decomposition_create(mrc_mat_comm(mat));

  sub->A = mrc_mat_create(MPI_COMM_SELF);
  mrc_mat_set_type(sub->A, "csr");

  sub->B = mrc_mat_create(MPI_COMM_SELF);
  mrc_mat_set_type(sub->B, "csr");
  
  sub->is_assembled = false;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_setup

static void
mrc_mat_csr_mpi_setup(struct mrc_mat *mat)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  sub->dc_row->n = mat->m;
  sub->dc_col->n = mat->n;

  mrc_decomposition_setup(sub->dc_row);
  mrc_decomposition_setup(sub->dc_col);

  // A is the diagonal block, so use local sizes
  mrc_mat_set_param_int(sub->A, "m", sub->dc_row->n);
  mrc_mat_set_param_int(sub->A, "n", sub->dc_col->n);
  mrc_mat_set_param_bool(sub->A, "verbose", sub->verbose);
  mrc_mat_set_param_int(sub->A, "nr_initial_cols", sub->nr_initial_cols);
  mrc_mat_set_param_float(sub->A, "growth_factor", sub->growth_factor);
  mrc_mat_setup(sub->A);

  // B is the off diagonal block, so # of rows = local # of rows,
  // but number of columns for now we set to the global number of columns
  // (even though local columns will never be inserted here, but
  // rather into A
  mrc_mat_set_param_int(sub->B, "m", sub->dc_row->n);
  mrc_mat_set_param_int(sub->B, "n", sub->dc_col->N);
  mrc_mat_set_param_bool(sub->B, "verbose", sub->verbose);
  mrc_mat_set_param_int(sub->B, "nr_initial_cols", sub->nr_initial_cols);
  mrc_mat_set_param_float(sub->B, "growth_factor", sub->growth_factor);
  mrc_mat_setup(sub->B);

  mrc_mat_setup_super(mat);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_destroy

static void
mrc_mat_csr_mpi_destroy(struct mrc_mat *mat)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

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
// mrc_mat_csr_mpi_add_value
//
// WARNING, all elements for any given row must be added contiguously!

static void
mrc_mat_csr_mpi_add_value(struct mrc_mat *mat, int row_idx, int col_idx, double val)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

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
// _mcsr_mpi_dump_mat

static void _mrc_unused
_csr_mpi_dump_mat(FILE *f, struct mrc_mat *mat_mpi, int which, bool ignore_identity)
{
  struct mrc_mat_csr_mpi *sub_mpi = mrc_mat_csr_mpi(mat_mpi);
  struct mrc_mat_csr *sub = NULL;
  int col_idx_off = 0;
  int m_off = sub_mpi->dc_col->off;
  
  if (which == 0) {
    sub = mrc_mat_csr(sub_mpi->A);
    col_idx_off = m_off;
  } else if (which == 1) {
    sub = mrc_mat_csr(sub_mpi->B);
  } else {
    assert(0);
  }

  mrc_fld_data_t *vals = mrc_vec_get_array(sub->vals);
  int *rows = mrc_vec_get_array(sub->rows);
  int *cols = mrc_vec_get_array(sub->cols);

  for (int row = 0; row < sub->nr_rows; row++) {
    for (int i=rows[row]; i < rows[row + 1]; i++) {
      int grow_idx = row + m_off;
      int gcol_idx = cols[i] + col_idx_off;
      mrc_fld_data_t val = vals[i];
      
      // don't print identity portion
      if (ignore_identity && grow_idx == gcol_idx && val == 1.0){
        continue;
      }
      fprintf(f, "%d\t\t%d\t\t%lg\n", grow_idx, gcol_idx, val);
    }
  }  
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_assemble

static void
mrc_mat_csr_mpi_assemble(struct mrc_mat *mat)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);
  int recv_len_max = 0;

  assert(!sub->is_assembled);
  
  mrc_mat_assemble(sub->A);
  mrc_mat_assemble(sub->B);

  int rank, size;
  MPI_Comm_rank(mrc_mat_comm(mat), &rank);
  MPI_Comm_size(mrc_mat_comm(mat), &size);

  // ////// DEBUGGING PRINT MATRICES TO FILE //////
  // static int assemble_count = 0;
  // char fname[80];
  // FILE *f;
  //  
  // sprintf(fname, "MAT_%s%d_A_%02d.txt", mat->obj.name, assemble_count, rank);
  // f = fopen(fname, "w");
  // mprintf(">> printing matrix %s @ assembly\n", fname);
  // _csr_mpi_dump_mat(f, mat, 0, false);
  // fclose(f);
  // 
  // sprintf(fname, "MAT_%s%d_B_%02d.txt", mat->obj.name, assemble_count, rank);
  // f = fopen(fname, "w");
  // mprintf(">> printing matrix %s @ assembly\n", fname);
  // _csr_mpi_dump_mat(f, mat, 1, false);
  // fclose(f);
  // assemble_count++;
  // ////// DEBUGGING PRINT MATRICES TO FILE //////  

  // create full ownership info on each proc
  int *col_offs_by_rank = calloc(size, sizeof(col_offs_by_rank));
  MPI_Allgather(&mat->n, 1, MPI_INT, col_offs_by_rank, 1, MPI_INT,
                mrc_mat_comm(mat));
  for (int r = 1; r < size; r++) {
    col_offs_by_rank[r] += col_offs_by_rank[r-1];
  }

  // set up map that maps non-zero column indices in B to a compacted index
  int N_cols = sub->dc_col->N;
  int *col_map = malloc(N_cols * sizeof(*col_map));
  for (int col = 0; col < N_cols; col++) {
    col_map[col] = -1;
  }

  // map each column with non-zero entries to an unique (per rank) entry in
  // the compact counting, starting at zero
  int *col_map_cnt_by_rank = calloc(size, sizeof(*col_map_cnt_by_rank));
  int col_map_cnt = 0;
  
  // FIXME: csr specific
  struct mrc_mat_csr *sub_B = mrc_mat_csr(sub->B);
  int *b_rows = mrc_vec_get_array(sub_B->rows);
  int *b_cols = mrc_vec_get_array(sub_B->cols);
  for (int row=0; row < sub_B->nr_rows; row++) {
    for (int i=b_rows[row]; i < b_rows[row + 1]; i++) {
      int col = b_cols[i];
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
  // FIXME: csr specific
  for (int row=0; row < sub_B->nr_rows; row++) {
    for (int i=b_rows[row]; i < b_rows[row + 1]; i++) {
      int col_idx = b_cols[i];
      // mprintf("> resetting B column: row %d  col %d -> %d\n",
      //         row, b_cols[i], col_map[col_idx]);
      b_cols[i] = col_map[col_idx];
    }
  }  
  sub->B->n = col_map_cnt;  // FIXME: this seems hacky, but doesn't break anything?
  free(col_map);
  mrc_vec_put_array(sub_B->rows, b_rows);
  mrc_vec_put_array(sub_B->cols, b_cols);
  b_rows = NULL;
  b_cols = NULL;  

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
      if(sub->recv_len[n] > recv_len_max) {
        recv_len_max = sub->recv_len[n];
      }
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
  // mprintf("%p alloc sub->req %p %d %d\n",
  //         mat, sub->req, sub->n_sends, sub->n_recvs);
  for (int n = 0; n < sub->n_sends; n++) {
    MPI_Irecv(&sub->send_len[n], 1, MPI_INT, MPI_ANY_SOURCE, 0, mrc_mat_comm(mat),
              &sub->req[n]);
  }

  for (int n = 0; n < sub->n_recvs; n++) {
    MPI_Isend(&sub->recv_len[n], 1, MPI_INT, sub->recv_src[n], 0, mrc_mat_comm(mat),
              &sub->req[sub->n_sends + n]);
  }

  MPI_Status *status = calloc(sub->n_sends + sub->n_recvs, sizeof(*status));
  MPI_Waitall(sub->n_sends + sub->n_recvs, sub->req, status);
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
  MPI_Waitall(sub->n_sends + sub->n_recvs, sub->req, MPI_STATUSES_IGNORE);

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
  
  // compile some rudimentary stats on what the matrix looks like
  if (sub->verbose) {
    int size;
    MPI_Comm_size(mrc_mat_comm(mat), &size);
    int commu_size_max = 0, recv_len_avg = 0;
    
    if (sub->n_recvs > 0) {
      recv_len_avg = col_map_cnt / sub->n_recvs;
    } else {
      recv_len_avg = 0;
    }
    mprintf("Sparse Mat Apply needs msgs from %d procs of size: %d (max),"
            " %d (avg), %d (total)\n", sub->n_recvs, recv_len_max,
            recv_len_avg, col_map_cnt);
    MPI_Reduce(&recv_len_max, &commu_size_max, 1, MPI_INT, MPI_MAX, 0,
               mrc_mat_comm(mat));
    mpi_printf(mrc_mat_comm(mat),
               "Sparse Mat max proc->proc communication: %d values\n",
               commu_size_max);

    // now print some stats on % of values that are non-local
    double pct_non_local = 0.0;  // % of matrix that's non-local
    double pct_non_local_max = 0.0, pct_non_local_sum = 0.0;
    double total_vals = (mrc_mat_csr(sub->A)->nr_vals +
                         mrc_mat_csr(sub->B)->nr_vals);
    if (total_vals > 0.0) {
      pct_non_local = ((double)mrc_mat_csr(sub->B)->nr_vals / total_vals);
    } else {
      pct_non_local = 0.0;
    }
    // mprintf("percent non-local:: %lg\n", 100.0 * pct_non_local);
    MPI_Reduce(&pct_non_local, &pct_non_local_max, 1, MPI_DOUBLE, MPI_MAX, 0,
               mrc_mat_comm(mat));
    MPI_Reduce(&pct_non_local, &pct_non_local_sum, 1, MPI_DOUBLE, MPI_SUM, 0,
               mrc_mat_comm(mat));
    mpi_printf(mrc_mat_comm(mat),
               "Sparse Mat percent non-local: %lg (max)  %lg (avg)\n",
               100.0 * pct_non_local_max, 100.0 * pct_non_local_sum / size);
  }
  
  sub->is_assembled = true;
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_gather_xc_start
//
// x can change before gather_xc_finish

static void
mrc_mat_csr_mpi_gather_xc_start(struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  assert(sub->is_assembled);
  assert(sub->_recv_vec == NULL && sub->_recv_buf == NULL);
  
  int *send_map_p;
  mrc_fld_data_t *pp;
  mrc_fld_data_t *x_arr = mrc_vec_get_array(x);
  
  sub->_recv_vec = sub->x_nl;
  sub->x_nl = NULL;  // you shouldn't touch x_nl... for real ;-)
  sub->_recv_buf = mrc_vec_get_array(sub->_recv_vec);
  
  // actual communication
  pp = sub->_recv_buf;
  for (int n = 0; n < sub->n_recvs; n++) {
    MPI_Irecv(pp, sub->recv_len[n], MPI_MRC_FLD_DATA_T, sub->recv_src[n], 2,
              mrc_mat_comm(mat), &sub->req[n]);
    pp += sub->recv_len[n];
  }
  assert(pp - sub->_recv_buf == mrc_vec_len(sub->_recv_vec));

  send_map_p = sub->send_map;
  pp = sub->send_buf;
  for (int n = 0; n < sub->n_sends; n++) {
    for (int i = 0; i < sub->send_len[n]; i++) {
      pp[i] = x_arr[send_map_p[i]];
    }
    MPI_Isend(pp, sub->send_len[n], MPI_MRC_FLD_DATA_T, sub->send_dst[n], 2,
              mrc_mat_comm(mat), &sub->req[sub->n_recvs + n]);
    send_map_p += sub->send_len[n];
    pp += sub->send_len[n];
  }

  mrc_vec_put_array(x, x_arr);
  // sub->_recv_buf is put on gather_xc_finish
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_gather_xc_finish

static void
mrc_mat_csr_mpi_gather_xc_finish(struct mrc_mat *mat)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);
  assert(sub->_recv_vec != NULL && sub->_recv_buf != NULL);
  
  MPI_Waitall(sub->n_sends + sub->n_recvs, sub->req, MPI_STATUSES_IGNORE);
  
  mrc_vec_put_array(sub->_recv_vec, sub->_recv_buf);
  sub->_recv_buf = NULL;
  sub->x_nl = sub->_recv_vec;
  sub->_recv_vec = NULL;  
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_apply
// y = mat * x

static void
mrc_mat_csr_mpi_apply(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  mrc_mat_csr_mpi_gather_xc_start(mat, x);
  mrc_mat_apply(y, sub->A, x);
  mrc_mat_csr_mpi_gather_xc_finish(mat);
  mrc_mat_apply_add(y, sub->B, sub->x_nl);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_apply_in_place
// x = mat * x

static void
mrc_mat_csr_mpi_apply_in_place(struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  // make a profiler for each unique matrix
#ifndef _NR_CSR_MPI_AIP_PROFS
#define _NR_CSR_MPI_AIP_PROFS (10)
#endif
  static void *mats[_NR_CSR_MPI_AIP_PROFS] = { NULL };
  static int apply_profs[_NR_CSR_MPI_AIP_PROFS] = { 0 };
  static int commu_profs[_NR_CSR_MPI_AIP_PROFS] = { 0 };
  static char name_apply[_NR_CSR_MPI_AIP_PROFS][20];
  static char name_commu[_NR_CSR_MPI_AIP_PROFS][20];
  static int pr_apply = 0;
  static int pr_commu = 0;
  static int pr_all = 0;
  static int pr_all_commu = 0;
  
  if (sub->do_profiling) {
    if (!pr_all) {
      pr_all = prof_register("csr_apply_in_place", 0, 0, 0);
      pr_all_commu = prof_register("csr_apply_in_place_commu", 0, 0, 0);
    }
    
    for(int i=0; i < _NR_CSR_MPI_AIP_PROFS; i++) {
      if (mats[i] == NULL) {
        // make new profilers
        snprintf(name_apply[i], 20, "csr_aip%d", i);
        snprintf(name_commu[i], 20, "csr_aip_commu%d", i);
        // mprintf("creating profiler: %s\n", name_apply[i]);
        // mprintf("creating profiler: %s\n", name_commu[i]);
        pr_apply = prof_register(name_apply[i], 0, 0, 0);
        pr_commu = prof_register(name_commu[i], 0, 0, 0);
        mats[i] = (void*)mat;
        apply_profs[i] = pr_apply;
        commu_profs[i] = pr_commu;
        break;
      } else if ((void*)mat == mats[i]) {
        // use existing profiler
        pr_apply = apply_profs[i];
        pr_commu = commu_profs[i];
        break;
      }
    }
    if (pr_apply == 0) mprintf("Note: Already profiling %d matrix apply_in_place's; "
                               "not profiling this one.\n", _NR_CSR_MPI_AIP_PROFS);
  }
  
  if (pr_all > 0) prof_start(pr_all);
  if (pr_apply > 0) prof_start(pr_apply);
  mrc_mat_csr_mpi_gather_xc_start(mat, x);
  mrc_mat_apply(x, sub->A, x);
  
  if (pr_all_commu > 0) prof_start(pr_all_commu);
  if (pr_commu > 0) prof_start(pr_commu);
  mrc_mat_csr_mpi_gather_xc_finish(mat);
  if (pr_commu > 0) prof_stop(pr_commu);
  if (pr_all_commu > 0) prof_stop(pr_all_commu);
  
  mrc_mat_apply_add(x, sub->B, sub->x_nl);
  
  if (pr_apply > 0) prof_stop(pr_apply);
  if (pr_all > 0) prof_stop(pr_all);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_apply_add
// y = mat * x + y

static void
mrc_mat_csr_mpi_apply_add(struct mrc_vec *y, struct mrc_mat *mat, struct mrc_vec *x)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  mrc_mat_csr_mpi_gather_xc_start(mat, x);
  mrc_mat_apply_add(y, sub->A, x);
  mrc_mat_csr_mpi_gather_xc_finish(mat);
  mrc_mat_apply_add(y, sub->B, sub->x_nl);
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_apply_general
// z = alpha * mat * x + beta * y

static void
mrc_mat_csr_mpi_apply_general(struct mrc_vec *z, double alpha,
                              struct mrc_mat *mat, struct mrc_vec *x,
                              double beta, struct mrc_vec *y)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  mrc_mat_csr_mpi_gather_xc_start(mat, x);
  mrc_mat_apply_general(z, alpha, sub->A, x, beta, y);
  mrc_mat_csr_mpi_gather_xc_finish(mat);
  mrc_mat_apply_general(z, alpha, sub->B, sub->x_nl, 1.0, z);
}


// ----------------------------------------------------------------------
// mrc_mat_csr_mpi_print

static void
mrc_mat_csr_mpi_print(struct mrc_mat *mat)
{
  struct mrc_mat_csr_mpi *sub = mrc_mat_csr_mpi(mat);

  mprintf("csr_mpi sub-matrix A:\n");
  mrc_mat_print(sub->A);
  mprintf("\n");

  mprintf("csr_mpi sub-matrix B:\n");
  mrc_mat_print(sub->B);
  mprintf("\n");
}

// ----------------------------------------------------------------------
// mrc_mat_csr_mpi description
#define VAR(x) (void *)offsetof(struct mrc_mat_csr_mpi, x)
static struct param mrc_mat_csr_mpi_descr[] = {
  { "verbose"           , VAR(verbose)           , PARAM_BOOL(false)    },
  { "do_profiling"      , VAR(do_profiling)      , PARAM_BOOL(false)    },
  { "nr_initial_cols"   , VAR(nr_initial_cols)   , PARAM_INT(1),
  .help = "How many empty columnts to use for each new row" },
  { "growth_factor"     , VAR(growth_factor)     , PARAM_FLOAT(0.5),
  .help = "Fraction of current nr_cols to add to rows" },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_mat subclass "csr_mpi"

struct mrc_mat_ops mrc_mat_csr_mpi_ops = {
  .name                  = "csr_mpi",
  .size                  = sizeof(struct mrc_mat_csr_mpi),
  .param_descr           = mrc_mat_csr_mpi_descr,
  .create                = mrc_mat_csr_mpi_create,
  .setup                 = mrc_mat_csr_mpi_setup,
  .destroy               = mrc_mat_csr_mpi_destroy,
  .add_value             = mrc_mat_csr_mpi_add_value,
  .assemble              = mrc_mat_csr_mpi_assemble,
  .apply                 = mrc_mat_csr_mpi_apply,
  .apply_in_place        = mrc_mat_csr_mpi_apply_in_place,
  .apply_add             = mrc_mat_csr_mpi_apply_add,
  .apply_general         = mrc_mat_csr_mpi_apply_general,
  .print                 = mrc_mat_csr_mpi_print,
};
