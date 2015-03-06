
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

  int *col_map;
  int col_map_cnt;
  int *rev_col_map;

  struct mrc_fld *xc;
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

  free(sub->col_map);
  free(sub->rev_col_map);

  mrc_fld_destroy(sub->xc);
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

  // set up map that maps non-zero column indices in B to a compacted index
  struct mrc_mat_mcsr *sub_B = mrc_mat_mcsr(sub->B);

  sub->col_map = malloc(sub->N * sizeof(*sub->col_map));
  int *col_map = sub->col_map;

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
  sub->col_map_cnt = col_map_cnt;

  for (int col = 0; col < sub->N; col++) {
    mprintf("map %d -> %d\n", col, col_map[col]);
  }

  // set up reverse map, mapping compacted indices back to original
  // global column indices
  sub->rev_col_map = malloc(sub->col_map_cnt * sizeof(*sub->rev_col_map));
  int *rev_col_map = sub->rev_col_map;

  for (int col = 0; col < sub->N; col++) {
    rev_col_map[sub->col_map[col]] = col;
  }
  for (int i = 0; i < sub->col_map_cnt; i++) {
    mprintf("rev map %d -> %d\n", i, rev_col_map[i]);
  }

  // update column indices in B to refer to compacted indices
  for (int row = 0; row < sub_B->nr_rows; row++) {
    for (int entry = sub_B->rows[row].first_entry;
	 entry < sub_B->rows[row + 1].first_entry; entry++) {
      int col_idx = sub_B->entries[entry].idx;
      sub_B->entries[entry].idx = sub->col_map[col_idx];
    }
  }

  // xc field compacted non-local part of right hand side vector x
  sub->xc = mrc_fld_create(mrc_mat_comm(mat));
  mrc_fld_set_type(sub->xc, FLD_TYPE);
  mrc_fld_set_param_int_array(sub->xc, "dims", 1, (int[1]) { sub->col_map_cnt });
  mrc_fld_setup(sub->xc);
}

static void
mrc_mat_mcsr_mpi_gather_xc(struct mrc_mat *mat, struct mrc_fld *x, struct mrc_fld *xc)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  // xg field -- FIXME, will be unneeded when we do proper communication
  struct mrc_fld *xg = mrc_fld_create(mrc_mat_comm(mat));
  mrc_fld_set_type(xg, FLD_TYPE);
  mrc_fld_set_param_int_array(xg, "dims", 1, (int[1]) { sub->M });
  mrc_fld_setup(xg);

  int rank, size;
  MPI_Comm_rank(mrc_mat_comm(mat), &rank);
  MPI_Comm_size(mrc_mat_comm(mat), &size);

  // create full ownership info on each proc
  int *col_offs_by_rank = calloc(size, sizeof(*col_offs_by_rank));
  MPI_Allgather(&mat->n, 1, MPI_INT, col_offs_by_rank, 1, MPI_INT,
		mrc_mat_comm(mat));
  for (int r = 1; r < size; r++) {
    col_offs_by_rank[r] += col_offs_by_rank[r-1];
  }
  /* for (int r = 0; r < size; r++) { */
  /*   mprintf("col_offs_by_rank[%d] = %d\n", r, col_offs_by_rank[r]); */
  /* } */

  // for each rank, find how many columns we need to receive from that rank
  int n_recvs = 0;
  int *recv_cnt_by_rank = calloc(size, sizeof(*recv_cnt_by_rank));
  int r = 0;
  for (int i = 0; i < sub->col_map_cnt; i++) {
    int col_idx = sub->rev_col_map[i];
    // if this col_idx isn't in the same proc's range as last time,
    // start search over
    if (r > 0 && col_idx < col_offs_by_rank[r-1]) {
      r = 0;
    }
    for (; r < size; r++) {
      if (col_idx < col_offs_by_rank[r]) {
	if (recv_cnt_by_rank[r]++ == 0) {
	  n_recvs++;
	} 
	break;
      }
    }
    assert(r < size);
  }
  /* mprintf("n_recvs = %d\n", n_recvs); */
  /* for (int r = 0; r < size; r++) { */
  /*   mprintf("recv_cnt_by_rank[%d] = %d\n", r, recv_cnt_by_rank[r]); */
  /* } */

  // we shouldn't have any columns to be received from our own rank,
  // because those are put into sub->A in the first place
  assert(recv_cnt_by_rank[rank] == 0);

  // compact recv_cnt_by_rank[]
  int *recv_len = calloc(n_recvs, sizeof(*recv_len));
  int *recv_src = calloc(n_recvs, sizeof(*recv_src));
  n_recvs = 0;
  for (int r = 0; r < size; r++) {
    if (recv_cnt_by_rank[r] > 0) {
      recv_len[n_recvs] = recv_cnt_by_rank[r];
      recv_src[n_recvs] = r;
      n_recvs++;
    }
  }
  free(recv_cnt_by_rank);

  // find out how many procs this proc needs to send messages to
  // (by finding this info for all procs first, then picking the
  // current rank's value)
  int *recv_flg_by_rank = calloc(size, sizeof(*recv_flg_by_rank));
  for (int n = 0; n < n_recvs; n++) {
    recv_flg_by_rank[recv_src[n]] = 1;
  }
  
  int *send_flg_by_rank = calloc(size, sizeof(*send_flg_by_rank));
  MPI_Allreduce(recv_flg_by_rank, send_flg_by_rank, size, MPI_INT, MPI_SUM,
		mrc_mat_comm(mat));
  int n_sends = send_flg_by_rank[rank];
  
  free(recv_flg_by_rank);
  free(send_flg_by_rank);

  // find compacted info on sends from the local proc,
  // by communication from those procs that expect to receive
  // those sends
  int *send_len = calloc(n_sends, sizeof(*send_len));
  int *send_dst = calloc(n_sends, sizeof(*send_dst));
  MPI_Request *req = calloc(n_sends + n_recvs, sizeof(*req));
  MPI_Status *status = calloc(n_sends + n_recvs, sizeof(*status));
  for (int n = 0; n < n_sends; n++) {
    MPI_Irecv(&send_len[n], 1, MPI_INT, MPI_ANY_SOURCE, 0, mrc_mat_comm(mat),
	      &req[n]);
  }

  for (int n = 0; n < n_recvs; n++) {
    MPI_Isend(&recv_len[n], 1, MPI_INT, recv_src[n], 0, mrc_mat_comm(mat),
	      &req[n_sends + n]);
  }

  MPI_Waitall(n_sends + n_recvs, req, status);
  for (int n = 0; n < n_sends; n++) {
    send_dst[n] = status[n].MPI_SOURCE;
  }

  for (int n = 0; n < n_recvs; n++) {
    mprintf("recv_len[%d] = %d src = %d\n", n, recv_len[n], recv_src[n]);
  }
  for (int n = 0; n < n_sends; n++) {
    mprintf("send_len[%d] = %d dst = %d\n", n, send_len[n], send_dst[n]);
  }

  // prepare actual send / recv buffers
  int send_buf_size = 0;
  for (int n = 0; n < n_sends; n++) {
    send_buf_size += send_len[n];
  }

  int recv_buf_size = 0;
  for (int n = 0; n < n_recvs; n++) {
    recv_buf_size += recv_len[n];
  }
  assert(recv_buf_size == sub->col_map_cnt);

  // now create a map on how to fill the send buffers
  int *send_map = calloc(send_buf_size, sizeof(*send_map));
  int *p = send_map;
  for (int n = 0; n < n_sends; n++) {
    MPI_Irecv(p, send_len[n], MPI_INT, send_dst[n], 1, mrc_mat_comm(mat),
	      &req[n]);
    p += send_len[n];
  }
  p = sub->rev_col_map;
  for (int n = 0; n < n_recvs; n++) {
    MPI_Isend(p, recv_len[n], MPI_INT, recv_src[n], 1, mrc_mat_comm(mat),
	      &req[n_sends + n]);
    p += recv_len[n];
  }

  MPI_Waitall(n_sends + n_recvs, req, MPI_STATUSES_IGNORE);

  p = send_map;
  for (int n = 0; n < n_sends; n++) {
    for (int i = 0; i < send_len[n]; i++) {
      if (rank > 0) {
	p[i] -= col_offs_by_rank[rank - 1];
      }
      assert(p[i] >= 0 && p[i] < mat->n);
      mprintf("map send %d: [%d] <- %d\n", n, i, p[i]);
    }
    p += send_len[n];
  }
  free(col_offs_by_rank);

  // actual communication
  mrc_fld_data_t *send_buf = calloc(send_buf_size, sizeof(*send_buf));
  mrc_fld_data_t *recv_buf = xc->_arr;

  mrc_fld_data_t *pp = send_buf;
  for (int n = 0; n < n_sends; n++) {
    for (int i = 0; i < send_len[n]; i++) {
      pp[i] = MRC_D1(x, send_map[i]);
    }
    MPI_Isend(pp, send_len[n], MPI_MRC_FLD_DATA_T, send_dst[n], 2, mrc_mat_comm(mat),
	      &req[n]);
    pp += send_len[n];
  }
  pp = recv_buf;
  for (int n = 0; n < n_recvs; n++) {
    MPI_Irecv(pp, recv_len[n], MPI_MRC_FLD_DATA_T, recv_src[n], 2, mrc_mat_comm(mat),
	      &req[n]);
    pp += recv_len[n];
  }
  MPI_Waitall(n_sends + n_recvs, req, MPI_STATUSES_IGNORE);

  free(req);
  free(status);

  // old
  MPI_Allgather(x->_arr, x->_len, MPI_MRC_FLD_DATA_T,
		xg->_arr, x->_len, MPI_MRC_FLD_DATA_T, mrc_mat_comm(mat));

  for (int i = 0; i < sub->col_map_cnt; i++) {
    mprintf("!!! %g %g\n", MRC_D1(xc, i), MRC_D1(xg, sub->rev_col_map[i]));
    MRC_D1(xc, i) = MRC_D1(xg, sub->rev_col_map[i]);
    mprintf("xc[%d] = %g\n", i, MRC_D1(xc, i));
  }

  mrc_fld_destroy(xg);
}

// ----------------------------------------------------------------------
// mrc_mat_mcsr_mpi_apply

static void
mrc_mat_mcsr_mpi_apply(struct mrc_fld *y, struct mrc_mat *mat, struct mrc_fld *x)
{
  struct mrc_mat_mcsr_mpi *sub = mrc_mat_mcsr_mpi(mat);

  mrc_mat_apply(y, sub->A, x);
  mrc_mat_mcsr_mpi_gather_xc(mat, x, sub->xc);
  mrc_mat_apply_add(y, sub->B, sub->xc);
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

