
#include <mrc_io_private.h>
#include <mrc_params.h>
#include "mrc_io_xdmf_lib.h"
#include <mrc_redist.h>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// ======================================================================

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

struct xdmf_file {
  hid_t h5_file;
  list_t xdmf_spatial_list;
};

struct xdmf {
  int slab_dims[3];
  int slab_off[3];
  struct xdmf_file file;
  struct xdmf_temporal *xdmf_temporal;
  bool use_independent_io;
  char *romio_cb_write;
  char *romio_ds_write;
  int nr_writers;
  MPI_Comm comm_writers; //< communicator for only the writers
  int *writers;          //< rank (in mrc_io comm) for each writer
  int is_writer;         //< this rank is a writer
  char *mode;            //< open mode, "r" or "w"
};

#define VAR(x) (void *)offsetof(struct xdmf, x)
static struct param xdmf_collective_descr[] = {
  { "use_independent_io"     , VAR(use_independent_io)      , PARAM_BOOL(false)      },
  { "nr_writers"             , VAR(nr_writers)              , PARAM_INT(1)           },
  { "romio_cb_write"         , VAR(romio_cb_write)          , PARAM_STRING(NULL)     },
  { "romio_ds_write"         , VAR(romio_ds_write)          , PARAM_STRING(NULL)     },
  { "slab_dims"              , VAR(slab_dims)               , PARAM_INT3(0, 0, 0)    },
  { "slab_off"               , VAR(slab_off)                , PARAM_INT3(0, 0, 0)    },
  {},
};
#undef VAR

#define to_xdmf(io) mrc_to_subobj(io, struct xdmf)

// ----------------------------------------------------------------------
// xdmf_collective_setup

static void
xdmf_collective_setup(struct mrc_io *io)
{
  mrc_io_setup_super(io);

  struct xdmf *xdmf = to_xdmf(io);

  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 7];
  sprintf(filename, "%s/%s.xdmf", io->par.outdir, io->par.basename);
  xdmf->xdmf_temporal = xdmf_temporal_create(filename);

#ifndef H5_HAVE_PARALLEL
  assert(xdmf->nr_writers == 1);
#endif
  
  if (xdmf->nr_writers > io->size) {
    xdmf->nr_writers = io->size;
  }
  xdmf->writers = calloc(xdmf->nr_writers, sizeof(*xdmf->writers));
  // setup writers, just use first nr_writers ranks,
  // could do something fancier in the future
  for (int i = 0; i < xdmf->nr_writers; i++) {
    xdmf->writers[i] = i;
    if (i == io->rank)
      xdmf->is_writer = 1;
  }
  MPI_Comm_split(mrc_io_comm(io), xdmf->is_writer, io->rank, &xdmf->comm_writers);
}

// ----------------------------------------------------------------------
// xdmf_collective_destroy

static void
xdmf_collective_destroy(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);
  
  free(xdmf->writers);
  if (xdmf->comm_writers) {
    MPI_Comm_free(&xdmf->comm_writers);
  }

  if (xdmf->xdmf_temporal) {
    xdmf_temporal_destroy(xdmf->xdmf_temporal);
  }
}

// ----------------------------------------------------------------------
// xdmf_collective_open

static void
xdmf_collective_open(struct mrc_io *io, const char *mode)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  xdmf->mode = strdup(mode);
  //  assert(strcmp(mode, "w") == 0);

  char filename[strlen(io->par.outdir) + strlen(io->par.basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", io->par.outdir, io->par.basename,
	  io->step, 0);

  if (xdmf->is_writer) {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info;
    MPI_Info_create(&info);
    if (xdmf->romio_cb_write) {
      MPI_Info_set(info, "romio_cb_write", xdmf->romio_cb_write);
    }
    if (xdmf->romio_ds_write) {
      MPI_Info_set(info, "romio_ds_write", xdmf->romio_ds_write);
    }
#ifdef H5_HAVE_PARALLEL
    H5Pset_fapl_mpio(plist, xdmf->comm_writers, info);
#endif
    if (strcmp(mode, "w") == 0) {
      file->h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    } else if (strcmp(mode, "r") == 0) {
      file->h5_file = H5Fopen(filename, H5F_ACC_RDONLY, plist);
    } else {
      assert(0);
    }
    H5Pclose(plist);
    MPI_Info_free(&info);
  }
  xdmf_spatial_open(&file->xdmf_spatial_list);
}

// ----------------------------------------------------------------------
// xdmf_collective_close

static void
xdmf_collective_close(struct mrc_io *io)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;

  xdmf_spatial_close(&file->xdmf_spatial_list, io, xdmf->xdmf_temporal);
  if (xdmf->is_writer) {
    H5Fclose(file->h5_file);
    memset(file, 0, sizeof(*file));
  }
  free(xdmf->mode);
  xdmf->mode = NULL;
}

static void
xdmf_collective_write_attr(struct mrc_io *io, const char *path, int type,
		const char *name, union param_u *pv)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;
  
  if (!xdmf->is_writer) {
    // FIXME? should check whether the attribute is the same on every proc?
    return;
  }
  hid_t group;
  if (H5Lexists(file->h5_file, path, H5P_DEFAULT) > 0) {
    group = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group);
  } else {
    group = H5Gcreate(file->h5_file, path, H5P_DEFAULT, H5P_DEFAULT,
		      H5P_DEFAULT); H5_CHK(group);
  }

  switch (type) {
  case PT_SELECT:
  case PT_INT:
  case MRC_VAR_INT:
    ierr = H5LTset_attribute_int(group, ".", name, &pv->u_int, 1); CE;
    break;
  case PT_BOOL: 
  case MRC_VAR_BOOL: {
    int val = pv->u_bool;
    ierr = H5LTset_attribute_int(group, ".", name, &val, 1); CE;
    break;
  }
  case PT_FLOAT:
  case MRC_VAR_FLOAT:
    ierr = H5LTset_attribute_float(group, ".", name, &pv->u_float, 1); CE;
    break;
  case PT_DOUBLE:
  case MRC_VAR_DOUBLE:
    ierr = H5LTset_attribute_double(group, ".", name, &pv->u_double, 1); CE;
    break;
  case PT_STRING:
    if (pv->u_string) {
      ierr = H5LTset_attribute_string(group, ".", name, pv->u_string); CE;
    } else {
      ierr = H5LTset_attribute_string(group, ".", name, "(NULL)"); CE;
    }
    break;
  case PT_INT3:
    ierr = H5LTset_attribute_int(group, ".", name, pv->u_int3, 3); CE;
    break;
  case PT_FLOAT3:
    ierr = H5LTset_attribute_float(group, ".", name, pv->u_float3, 3); CE;
    break;
  case PT_DOUBLE3:
  case MRC_VAR_DOUBLE3:
    ierr = H5LTset_attribute_double(group, ".", name, pv->u_double3, 3); CE;
    break;
  case PT_INT_ARRAY: {
    hsize_t dims = pv->u_int_array.nr_vals;
    hid_t dataspace_id = H5Screate_simple(1, &dims, NULL); H5_CHK(dataspace_id);
    hid_t attr_id = H5Acreate(group, name, H5T_NATIVE_INT, dataspace_id,
			      H5P_DEFAULT, H5P_DEFAULT); H5_CHK(attr_id);
    if (dims > 0) {
      ierr = H5Awrite(attr_id, H5T_NATIVE_INT, pv->u_int_array.vals); CE;
    }
    ierr = H5Sclose(dataspace_id); CE;
    ierr = H5Aclose(attr_id); CE;
    break;
  }
  case PT_PTR:
    break;
  default:
    mprintf("mrc_io_xdmf_collective: not writing attr '%s' (type %d)\n",
	    name, type);
    assert(0);
  }
  ierr = H5Gclose(group); CE;
}

static void
xdmf_collective_read_attr(struct mrc_io *io, const char *path, int type,
			  const char *name, union param_u *pv)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  // read on I/O procs
  if (xdmf->is_writer) {
    hid_t group = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group);
    switch (type) {
    case PT_SELECT:
    case PT_INT:
    case MRC_VAR_INT:
      ierr = H5LTget_attribute_int(group, ".", name, &pv->u_int); CE;
      break;
    case PT_BOOL:
    case MRC_VAR_BOOL: ;
      int val;
      ierr = H5LTget_attribute_int(group, ".", name, &val); CE;
      pv->u_bool = val;
      break;
    case PT_FLOAT:
      ierr = H5LTget_attribute_float(group, ".", name, &pv->u_float); CE;
      break;
    case PT_DOUBLE:
      ierr = H5LTget_attribute_double(group, ".", name, &pv->u_double); CE;
      break;
    case PT_STRING: ;
      hsize_t dims;
      H5T_class_t class;
      size_t sz;
      ierr = H5LTget_attribute_info(group, ".", name, &dims, &class, &sz); CE;
      pv->u_string = malloc(sz);
      ierr = H5LTget_attribute_string(group, ".", name, (char *)pv->u_string); CE;
      if (strcmp(pv->u_string, "(NULL)") == 0) {
      	free((char *) pv->u_string);
      	pv->u_string = NULL;
      }
      break;
    case PT_INT3:
      ierr = H5LTget_attribute_int(group, ".", name, pv->u_int3); CE;
      break;
    case PT_FLOAT3:
      ierr = H5LTget_attribute_float(group, ".", name, pv->u_float3); CE;
      break;
    case PT_DOUBLE3:
      ierr = H5LTget_attribute_double(group, ".", name, pv->u_double3); CE;
      break;
    case PT_INT_ARRAY: {
      int attr = H5Aopen(group, name, H5P_DEFAULT); H5_CHK(attr);
      H5A_info_t ainfo;
      ierr = H5Aget_info(attr, &ainfo); CE;
      ierr = H5Aclose(attr); CE;
      pv->u_int_array.nr_vals = ainfo.data_size / sizeof(int);
      pv->u_int_array.vals = calloc(pv->u_int_array.nr_vals, sizeof(int));
      ierr = H5LTget_attribute_int(group, ".", name, pv->u_int_array.vals); CE;
      break;
    }
    case PT_PTR:
      break;
    default:
      mprintf("mrc_io_xdmf_collective: not reading attr '%s' (type %d)\n", name, type);
      assert(0);
      break;
    }
    ierr = H5Gclose(group); CE;
  }

  int root = xdmf->writers[0];
  MPI_Comm comm = mrc_io_comm(io);
  switch (type) {
  case PT_SELECT:
  case PT_INT:
  case MRC_VAR_INT:
    MPI_Bcast(&pv->u_int, 1, MPI_INT, root, comm);
    break;
  case PT_BOOL:
  case MRC_VAR_BOOL: ;
    int val = pv->u_int;
    MPI_Bcast(&val, 1, MPI_INT, root, comm);
    pv->u_int = val;
    break;
  case PT_FLOAT:
    MPI_Bcast(&pv->u_float, 1, MPI_FLOAT, root, comm);
    break;
  case PT_DOUBLE:
    MPI_Bcast(&pv->u_double, 1, MPI_DOUBLE, root, comm);
    break;
  case PT_STRING: ;
    int len;
    if (io->rank == root) {
      len = pv->u_string ? strlen(pv->u_string) : 0;
    }
    MPI_Bcast(&len, 1, MPI_INT, root, comm);

    if (!len) {
      pv->u_string = NULL;
      break;
    }

    if (io->rank != root) {
      pv->u_string = malloc(len + 1);
    }
    // FIXME, u_string type should not be const
    MPI_Bcast((char *) pv->u_string, len + 1, MPI_CHAR, root, comm);
    break;
  case PT_INT3:
    MPI_Bcast(pv->u_int3, 3, MPI_INT, root, comm);
    break;
  case PT_FLOAT3:
    MPI_Bcast(pv->u_float3, 3, MPI_FLOAT, root, comm);
    break;
  case PT_DOUBLE3:
    MPI_Bcast(pv->u_double3, 3, MPI_DOUBLE, root, comm);
    break;
  case PT_INT_ARRAY:
    MPI_Bcast(&pv->u_int_array.nr_vals, 1, MPI_INT, root, comm);
    if (io->rank != root) {
      free(pv->u_int_array.vals);
      pv->u_int_array.vals = calloc(pv->u_int_array.nr_vals, sizeof(int));
    }
    MPI_Bcast(pv->u_int_array.vals, pv->u_int_array.nr_vals, MPI_INT, root, comm);
    break;
  case PT_PTR:
    break;
  default:
    mprintf("mrc_io_xdmf_collective: attr '%s' (type %d)\n", name, type);
    assert(0);
    break;
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// collective_m1_write_f1
// does the actual write of the f1 to the file
// only called on writer procs

static void
collective_m1_write_f1(struct mrc_io *io, const char *path, struct mrc_ndarray *nd1,
		       int m, const char *name, hid_t group0)
{
  struct xdmf *xdmf = to_xdmf(io);
  int ierr;

  double io_scale;
  int rc = mrc_ndarray_get_var_double(nd1, "io_scale", &io_scale);
  if (rc == 0) {
    mrc_ndarray_scale(nd1, io_scale);
  }
  assert(mrc_ndarray_data_type(nd1) == MRC_NT_FLOAT);
  hid_t group_fld = H5Gcreate(group0, name, H5P_DEFAULT,
			      H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);
  ierr = H5LTset_attribute_int(group_fld, ".", "m", &m, 1); CE;
  
  hid_t group = H5Gcreate(group_fld, "p0", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
  int i0 = 0;
  ierr = H5LTset_attribute_int(group, ".", "global_patch", &i0, 1); CE;

  hsize_t fdims[1] = { mrc_ndarray_dims(nd1)[0] };
  hid_t filespace = H5Screate_simple(1, fdims, NULL); H5_CHK(filespace);
  hid_t dset = H5Dcreate(group, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl); // FIXME, consolidate
#ifdef H5_HAVE_PARALLEL
  if (xdmf->use_independent_io) {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT); CE;
  } else {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); CE;
  }
#endif
  hid_t memspace;
  if (io->rank == xdmf->writers[0]) {
    memspace = H5Screate_simple(1, fdims, NULL);
  } else {
    memspace = H5Screate(H5S_NULL);
    H5Sselect_none(memspace);
    H5Sselect_none(filespace);
  }
  ierr = H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, nd1->arr); CE;
  
  ierr = H5Dclose(dset); CE;
  ierr = H5Sclose(memspace); CE;
  ierr = H5Sclose(filespace); CE;
  ierr = H5Pclose(dxpl); CE;

  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group_fld); CE;
}

struct collective_m1_ctx {
  int nr_patches;
  int nr_global_patches;
  int dim;
  int sw;
  int gdims[3];
  int np[3];
  MPI_Request *send_reqs;
  int nr_send_reqs;
  MPI_Request *recv_reqs;
  int nr_recv_reqs;
  char comp_name[100];
};

static void
collective_m1_send_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			 struct mrc_fld *m1, int m)
{
  struct xdmf *xdmf = to_xdmf(io);
  int dim = ctx->dim;

  assert(mrc_fld_data_type(m1) == MRC_NT_FLOAT);
  ctx->send_reqs = calloc(ctx->nr_patches, sizeof(*ctx->send_reqs));
  ctx->nr_send_reqs = 0;

  for (int p = 0; p < ctx->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(m1->_domain, p, &info);
    bool skip = false;
    for (int d = 0; d < 3; d++) {
      if (d != dim && info.off[d] != 0) {
	skip = true;
      }
    }
    if (skip) {
      continue;
    }
    
    // FIXME, should use intersection, probably won't work if slab_dims are actually smaller
    int ib = 0;
    if (info.off[dim] == 0) { // FIXME, -> generic code
      ib = xdmf->slab_off[dim];
    }
    int ie = info.ldims[dim];
    if (info.off[dim] + info.ldims[dim] == ctx->gdims[dim]) {
      ie = xdmf->slab_off[dim] + xdmf->slab_dims[dim] - info.off[dim];
    }
    //mprintf("send to %d tag %d len %d\n", xdmf->writers[0], info.global_patch, ie - ib);
    assert(ib < ie);
    MPI_Isend(&MRC_M1(m1, m, ib, p), ie - ib, MPI_FLOAT,
	      xdmf->writers[0], info.global_patch, mrc_io_comm(io),
	      &ctx->send_reqs[ctx->nr_send_reqs++]);
  }
}

static void
collective_m1_send_end(struct mrc_io *io, struct collective_m1_ctx *ctx)
{
  MPI_Waitall(ctx->nr_send_reqs, ctx->send_reqs, MPI_STATUSES_IGNORE);
  free(ctx->send_reqs);
}

static void
collective_m1_recv_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			 struct mrc_domain *domain, struct mrc_ndarray *nd1)
{
  struct xdmf *xdmf = to_xdmf(io);

  if (io->rank != xdmf->writers[0])
    return;

  int dim = ctx->dim;

  ctx->recv_reqs = calloc(ctx->np[dim], sizeof(*ctx->recv_reqs));
  ctx->nr_recv_reqs = 0;

  for (int gp = 0; gp < ctx->nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);
    bool skip = false;
    for (int d = 0; d < 3; d++) {
      if (d != dim && info.off[d] != 0) {
	skip = true;
      }
    }
    if (skip) {
      continue;
    }
    
    int ib = info.off[dim];
    if (ib == 0) {
      ib = xdmf->slab_off[dim];
    }
    int ie = info.off[dim] + info.ldims[dim];
    if (ie == ctx->gdims[dim]) {
      ie = xdmf->slab_off[dim] + xdmf->slab_dims[dim];
    }
    //mprintf("recv from %d tag %d len %d\n", info.rank, gp, ie - ib);
    MPI_Irecv(&MRC_S1(nd1, ib), ie - ib, MPI_FLOAT, info.rank,
	      gp, mrc_io_comm(io), &ctx->recv_reqs[ctx->nr_recv_reqs++]);
  }
  assert(ctx->nr_recv_reqs == ctx->np[dim]);
}

static void
collective_m1_recv_end(struct mrc_io *io, struct collective_m1_ctx *ctx)
{
  struct xdmf *xdmf = to_xdmf(io);

  if (io->rank != xdmf->writers[0])
    return;

  MPI_Waitall(ctx->nr_recv_reqs, ctx->recv_reqs, MPI_STATUSES_IGNORE);
  free(ctx->recv_reqs);
}

// ----------------------------------------------------------------------
// xdmf_collective_write_m1

static void
xdmf_collective_write_m1(struct mrc_io *io, const char *path, struct mrc_fld *m1)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  struct collective_m1_ctx ctx;
  int nr_comps = mrc_fld_nr_comps(m1);
  mrc_fld_get_param_int(m1, "dim", &ctx.dim);
  ctx.sw = m1->_sw.vals[0];
  mrc_domain_get_global_dims(m1->_domain, ctx.gdims);
  mrc_domain_get_nr_global_patches(m1->_domain, &ctx.nr_global_patches);
  mrc_domain_get_nr_procs(m1->_domain, ctx.np);
  mrc_domain_get_patches(m1->_domain, &ctx.nr_patches);
  int dim = ctx.dim;
  int slab_off_save, slab_dims_save;
  slab_off_save = xdmf->slab_off[dim];
  slab_dims_save = xdmf->slab_dims[dim];
  // FIXME
  if (!xdmf->slab_dims[dim]) {
    xdmf->slab_dims[dim] = ctx.gdims[dim] + 2 * ctx.sw;
    xdmf->slab_off[dim] = -ctx.sw;
  }

  if (xdmf->is_writer) {
    // we're creating the nd1 on all writers, but only fill and actually write
    // it on writers[0]
    struct mrc_ndarray *nd1 = mrc_ndarray_create(MPI_COMM_SELF);
    mrc_ndarray_set_param_int_array(nd1, "dims", 1, (int [1]) { xdmf->slab_dims[dim] });
    mrc_ndarray_set_param_int_array(nd1, "offs", 1, (int [1]) { xdmf->slab_off[dim]  });
    mrc_ndarray_setup(nd1);

    double io_scale;
    int rc = mrc_fld_get_var_double(m1, "io_scale", &io_scale);
    if (rc == 0) {
      mrc_ndarray_dict_add_double(nd1, "io_scale", io_scale);
    }

    hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group0);
    for (int m = 0; m < nr_comps; m++) {
      collective_m1_recv_begin(io, &ctx, m1->_domain, nd1);
      collective_m1_send_begin(io, &ctx, m1, m);
      collective_m1_recv_end(io, &ctx);
      collective_m1_write_f1(io, path, nd1, 0, mrc_fld_comp_name(m1, m), group0);
      collective_m1_send_end(io, &ctx);
    }
    ierr = H5Gclose(group0); CE;

    mrc_ndarray_destroy(nd1);
  } else { // not writer
    for (int m = 0; m < nr_comps; m++) {
      collective_m1_send_begin(io, &ctx, m1, m);
      collective_m1_send_end(io, &ctx);
    }
  }
  xdmf->slab_dims[dim] = slab_dims_save;
  xdmf->slab_off[dim] = slab_off_save;
}

// ----------------------------------------------------------------------

static void
collective_m1_read_recv_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			      struct mrc_fld *m1, int m)
{
  struct xdmf *xdmf = to_xdmf(io);

  ctx->recv_reqs = calloc(1 + ctx->nr_patches, sizeof(*ctx->recv_reqs));
  ctx->nr_recv_reqs = 0;

  MPI_Irecv(ctx->comp_name, 100, MPI_CHAR, xdmf->writers[0], 0, mrc_io_comm(io),
	    &ctx->recv_reqs[ctx->nr_recv_reqs++]);

  for (int p = 0; p < ctx->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(m1->_domain, p, &info);
    int ib = -ctx->sw;
    int ie = info.ldims[ctx->dim] + ctx->sw;
    //	mprintf("recv to %d tag %d\n", xdmf->writers[0], info.global_patch);
    MPI_Irecv(&MRC_M1(m1, m, ib, p), ie - ib, MPI_FLOAT,
	      xdmf->writers[0], info.global_patch, mrc_io_comm(io),
	      &ctx->recv_reqs[ctx->nr_recv_reqs++]);
  }
}

static void
collective_m1_read_recv_end(struct mrc_io *io, struct collective_m1_ctx *ctx,
			    struct mrc_fld *m1, int m)
{
  MPI_Waitall(ctx->nr_recv_reqs, ctx->recv_reqs, MPI_STATUSES_IGNORE);
  free(ctx->recv_reqs);
  mrc_fld_set_comp_name(m1, m, ctx->comp_name);
}

static void
collective_m1_read_send_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			      struct mrc_domain *domain, struct mrc_fld *f1, int m)
{
  struct xdmf *xdmf = to_xdmf(io);

  if (io->rank != xdmf->writers[0])
    return;

  int dim = ctx->dim;
  ctx->send_reqs = calloc(io->size + ctx->nr_global_patches, sizeof(*ctx->send_reqs));
  ctx->nr_send_reqs = 0;
  const char *comp_name = mrc_fld_comp_name(f1, m);
  assert(comp_name && strlen(comp_name) < 99);
  for (int r = 0; r < io->size; r++) {
    MPI_Isend((char *) comp_name, strlen(comp_name) + 1, MPI_CHAR,
	      r, 0, mrc_io_comm(io), &ctx->send_reqs[ctx->nr_send_reqs++]);
  }
  for (int gp = 0; gp < ctx->nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);
    int ib = info.off[dim] - ctx->sw;
    int ie = info.off[dim] + info.ldims[dim] + ctx->sw;
    //  mprintf("send from %d tag %d\n", info.rank, gp);
    MPI_Isend(&MRC_F1(f1, m, ib), ie - ib, MPI_FLOAT,
	      info.rank, gp, mrc_io_comm(io), &ctx->send_reqs[ctx->nr_send_reqs++]);
  }
}

static void
collective_m1_read_send_end(struct mrc_io *io, struct collective_m1_ctx *ctx)
{
  struct xdmf *xdmf = to_xdmf(io);

  if (io->rank != xdmf->writers[0])
    return;

  MPI_Waitall(ctx->nr_send_reqs, ctx->send_reqs, MPI_STATUSES_IGNORE);
  free(ctx->send_reqs);
}

struct read_m1_cb_data {
  struct mrc_io *io;
  struct mrc_fld *gfld;
  hid_t filespace;
  hid_t memspace;
  hid_t dxpl;
};

static herr_t
read_m1_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{
  struct read_m1_cb_data *data = op_data;
  struct mrc_fld *gfld = data->gfld;
  int *ib = gfld->_ghost_offs;
  int ierr;

  hid_t group_fld = H5Gopen(g_id, name, H5P_DEFAULT); H5_CHK(group_fld);
  int m;
  ierr = H5LTget_attribute_int(group_fld, ".", "m", &m); CE;
  mrc_fld_set_comp_name(data->gfld, m, name);
  hid_t group = H5Gopen(group_fld, "p0", H5P_DEFAULT); H5_CHK(group);

  hid_t dset = H5Dopen(group, "1d", H5P_DEFAULT); H5_CHK(dset);
  ierr = H5Dread(dset, H5T_NATIVE_FLOAT, data->memspace, data->filespace,
		 data->dxpl, &MRC_F1(gfld, m, ib[0])); CE;
  ierr = H5Dclose(dset); CE;
  
  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group_fld); CE;

  return 0;
}

// ----------------------------------------------------------------------
// xdmf_collective_read_m1

static void
xdmf_collective_read_m1(struct mrc_io *io, const char *path, struct mrc_fld *m1)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  assert(mrc_fld_data_type(m1) == MRC_NT_FLOAT);
  struct collective_m1_ctx ctx;
  int gdims[3];
  int nr_comps = mrc_fld_nr_comps(m1);
  mrc_fld_get_param_int(m1, "dim", &ctx.dim);
  ctx.sw = m1->_sw.vals[0];
  mrc_domain_get_global_dims(m1->_domain, gdims);
  mrc_domain_get_nr_global_patches(m1->_domain, &ctx.nr_global_patches);
  mrc_domain_get_patches(m1->_domain, &ctx.nr_patches);

  if (xdmf->is_writer) {
    struct mrc_fld *f1 = mrc_fld_create(MPI_COMM_SELF);
    mrc_fld_set_param_int_array(f1, "dims", 2, (int [2]) { gdims[ctx.dim], nr_comps });
    mrc_fld_set_param_int_array(f1, "sw", 2, (int [2]) { ctx.sw, 0 });
    mrc_fld_setup(f1);

    hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group0);

    hsize_t hgdims[1] = { mrc_fld_ghost_dims(f1)[0] };

    hid_t filespace = H5Screate_simple(1, hgdims, NULL); H5_CHK(filespace);
    hid_t memspace;
    if (io->rank == xdmf->writers[0]) {
      memspace = H5Screate_simple(1, hgdims, NULL); H5_CHK(memspace);
    } else {
      memspace = H5Screate(H5S_NULL);
      H5Sselect_none(memspace);
      H5Sselect_none(filespace);
    }
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl);
#ifdef H5_HAVE_PARALLEL
    if (xdmf->use_independent_io) {
      H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
    } else {
      H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    }
#endif

    struct read_m1_cb_data cb_data = {
      .io        = io,
      .gfld      = f1,
      .memspace  = memspace,
      .filespace = filespace,
      .dxpl      = dxpl,
    };
    
    hsize_t idx = 0;
    H5Literate_by_name(group0, ".", H5_INDEX_NAME, H5_ITER_INC, &idx,
		       read_m1_cb, &cb_data, H5P_DEFAULT);

    ierr = H5Pclose(dxpl); CE;
    ierr = H5Sclose(memspace); CE;
    ierr = H5Sclose(filespace); CE;

    ierr = H5Gclose(group0); CE;

    for (int m = 0; m < nr_comps; m++) {
      collective_m1_read_recv_begin(io, &ctx, m1, m);
      collective_m1_read_send_begin(io, &ctx, m1->_domain, f1, m);
      collective_m1_read_recv_end(io, &ctx, m1, m);
      collective_m1_read_send_end(io, &ctx);
    }
  } else { // not writer
    for (int m = 0; m < nr_comps; m++) {
      collective_m1_read_recv_begin(io, &ctx, m1, m);
      collective_m1_read_recv_end(io, &ctx, m1, m);
    }
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// writer_write_fld
// does the actual write of the partial fld to the file
// only called on writer procs

static void
writer_write_fld(struct mrc_redist *redist, struct mrc_io *io,
		 const char *path, struct mrc_ndarray *nd, int m,
		 struct mrc_fld *m3, struct xdmf_spatial *xs, hid_t group0)
{
  int ierr;

  char default_name[100];
  const char *compname;

  // If the comps aren't named just name them by their component number
  if ( !(compname = mrc_fld_comp_name(m3, m)) ) {
    sprintf(default_name, "_UNSET_%d", m);
    compname = (const char *) default_name;
  }

  xdmf_spatial_save_fld_info(xs, strdup(compname), strdup(path), false, mrc_fld_data_type(m3));

  hid_t group_fld = H5Gcreate(group0, compname, H5P_DEFAULT,
			      H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);
  ierr = H5LTset_attribute_int(group_fld, ".", "m", &m, 1); CE;
  
  hid_t group = H5Gcreate(group_fld, "p0", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
  int i0 = 0;
  ierr = H5LTset_attribute_int(group, ".", "global_patch", &i0, 1); CE;

  hsize_t fdims[3] = { redist->slab_dims[2], redist->slab_dims[1], redist->slab_dims[0] };
  hid_t filespace = H5Screate_simple(3, fdims, NULL); H5_CHK(filespace);
  hid_t dtype;
  switch (mrc_ndarray_data_type(nd)) {
  case MRC_NT_FLOAT: dtype = H5T_NATIVE_FLOAT; break;
  case MRC_NT_DOUBLE: dtype = H5T_NATIVE_DOUBLE; break;
  case MRC_NT_INT: dtype = H5T_NATIVE_INT; break;
  default: assert(0);
  }

  hid_t dset = H5Dcreate(group, "3d", dtype, filespace, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl);
#ifdef H5_HAVE_PARALLEL
  struct xdmf *xdmf = to_xdmf(io);
  if (xdmf->use_independent_io) {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT); CE;
  } else {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); CE;
  }
#endif
  const int *im = mrc_ndarray_dims(nd), *ib = mrc_ndarray_offs(nd);
  hsize_t mdims[3] = { im[2], im[1], im[0] };
  hsize_t foff[3] = { ib[2] - redist->slab_offs[2],
		      ib[1] - redist->slab_offs[1],
		      ib[0] - redist->slab_offs[0] };
  hid_t memspace = H5Screate_simple(3, mdims, NULL);
  ierr = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL,
			     mdims, NULL); CE;

  ierr = H5Dwrite(dset, dtype, memspace, filespace, dxpl, nd->arr); CE;
  
  ierr = H5Dclose(dset); CE;
  ierr = H5Sclose(memspace); CE;
  ierr = H5Sclose(filespace); CE;
  ierr = H5Pclose(dxpl); CE;

  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group_fld); CE;
}

// ----------------------------------------------------------------------
// xdmf_collective_write_m3

static void
xdmf_collective_write_m3(struct mrc_io *io, const char *path, struct mrc_fld *m3)
{
  struct xdmf *xdmf = to_xdmf(io);

  struct mrc_redist redist[1];
  mrc_redist_init(redist, m3->_domain, xdmf->slab_off, xdmf->slab_dims,
		  xdmf->nr_writers);

  struct xdmf_file *file = &xdmf->file;
  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
					      mrc_domain_name(m3->_domain));
  if (!xs) {
    xs = xdmf_spatial_create_m3_parallel(&file->xdmf_spatial_list,
					 mrc_domain_name(m3->_domain),
					 m3->_domain,
					 redist->slab_offs, redist->slab_dims, io);
  }

  // If we have an aos field, we need to get it as soa for collection and writing
  struct mrc_fld *m3_soa = m3;
  if (m3->_aos) {
    switch (mrc_fld_data_type(m3)) {
      case MRC_NT_FLOAT: m3_soa = mrc_fld_get_as(m3, "float"); break;
      case MRC_NT_DOUBLE: m3_soa = mrc_fld_get_as(m3, "double"); break;
      case MRC_NT_INT: m3_soa = mrc_fld_get_as(m3, "int"); break;
      default: assert(0);
    }
  }

  hid_t group0 = 0;

  if (xdmf->is_writer) {
    if (H5Lexists(file->h5_file, path, H5P_DEFAULT) > 0) {
      group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
    } else {
      assert(0); // FIXME, can this happen?
      group0 = H5Gcreate(file->h5_file, path, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group0);
    }
    int nr_1 = 1;
    H5LTset_attribute_int(group0, ".", "nr_patches", &nr_1, 1);
  }

  struct mrc_ndarray *nd = mrc_redist_get_ndarray(redist, m3_soa);

  for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
    mrc_redist_run(redist, nd, m3_soa, m);

    if (redist->is_writer) {
      writer_write_fld(redist, io, path, nd, m, m3, xs, group0);
    }
  }

  mrc_redist_put_ndarray(redist, nd);
  
  if (xdmf->is_writer) {
    H5Gclose(group0);
  }

  if (m3->_aos) {
    mrc_fld_put_as(m3_soa, m3);
  }

  mrc_redist_destroy(redist);
}

// ======================================================================

// ----------------------------------------------------------------------
// collective helper context

struct collective_m3_entry {
  struct mrc_fld *fld;
  int ilo[3];
  int ihi[3];
  int patch;
  int global_patch; //< also used as tag
  int rank; //< of peer
};

struct mrc_redist_read {
  struct collective_m3_entry *blocks;
  MPI_Request *reqs;
  int n;
};

struct collective_m3_ctx {
  int gdims[3];
  int slab_dims[3], slab_off[3];
  int nr_patches, nr_global_patches;
  int slow_dim;
  int slow_indices_per_writer;
  int slow_indices_rmndr;

  struct mrc_redist_read read_send;
  struct mrc_redist_read read_recv;
};

static void
collective_m3_init(struct mrc_io *io, struct collective_m3_ctx *ctx,
		   struct mrc_domain *domain)
{
  struct xdmf *xdmf = to_xdmf(io);

  mrc_domain_get_global_dims(domain, ctx->gdims);
  mrc_domain_get_patches(domain, &ctx->nr_patches);
  mrc_domain_get_nr_global_patches(domain, &ctx->nr_global_patches);
  for (int d = 0; d < 3; d++) {
    if (xdmf->slab_dims[d]) {
      ctx->slab_dims[d] = xdmf->slab_dims[d];
    } else {
      ctx->slab_dims[d] = ctx->gdims[d];
    }
    ctx->slab_off[d] = xdmf->slab_off[d];
  }
  ctx->slow_dim = 2;
  while (ctx->gdims[ctx->slow_dim] == 1) {
    ctx->slow_dim--;
  }
  assert(ctx->slow_dim >= 0);
  int total_slow_indices = ctx->slab_dims[ctx->slow_dim];
  ctx->slow_indices_per_writer = total_slow_indices / xdmf->nr_writers;
  ctx->slow_indices_rmndr = total_slow_indices % xdmf->nr_writers;
}

static void
get_writer_off_dims(struct collective_m3_ctx *ctx, int writer,
		    int *writer_off, int *writer_dims)
{
  for (int d = 0; d < 3; d++) {
    writer_dims[d] = ctx->slab_dims[d];
    writer_off[d] = ctx->slab_off[d];
  }
  writer_dims[ctx->slow_dim] = ctx->slow_indices_per_writer + (writer < ctx->slow_indices_rmndr);
  if (writer < ctx->slow_indices_rmndr) {
    writer_off[ctx->slow_dim] += (ctx->slow_indices_per_writer + 1) * writer;
  } else {
    writer_off[ctx->slow_dim] += ctx->slow_indices_rmndr +
      ctx->slow_indices_per_writer * writer;
  }
}

static void
collective_m3_send_setup(struct mrc_io *io, struct collective_m3_ctx *ctx,
			 struct mrc_domain *domain, struct mrc_fld *gfld)
{
  struct mrc_redist_read *send = &ctx->read_send;

  send->n = 0;
  for (int gp = 0; gp < ctx->nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);
    
    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     mrc_fld_ghost_offs(gfld), mrc_fld_ghost_dims(gfld));
    if (has_intersection)
      send->n++;
  }

  send->reqs = calloc(send->n, sizeof(*send->reqs));
  send->blocks = calloc(send->n, sizeof(*send->blocks));

  for (int i = 0, gp = 0; i < send->n; i++) {
    struct collective_m3_entry *block = &send->blocks[i];
    struct mrc_patch_info info;
    bool has_intersection;
    do {
      mrc_domain_get_global_patch_info(domain, gp++, &info);
      has_intersection = find_intersection(block->ilo, block->ihi, info.off, info.ldims,
					   mrc_fld_ghost_offs(gfld), mrc_fld_ghost_dims(gfld));
    } while (!has_intersection);
    block->patch = info.global_patch;
    block->rank = info.rank;
  }
}

static void
collective_m3_send_begin(struct mrc_io *io, struct collective_m3_ctx *ctx,
			 struct mrc_domain *domain, struct mrc_fld *gfld)
{
  struct mrc_redist_read *send = &ctx->read_send;

  collective_m3_send_setup(io, ctx, domain, gfld);

  for (int i = 0; i < send->n; i++) {
    struct collective_m3_entry *block = &send->blocks[i];
    int *ilo = block->ilo, *ihi = block->ihi;
    struct mrc_fld *fld = mrc_fld_create(MPI_COMM_NULL);
    mrc_fld_set_type(fld, mrc_fld_type(gfld));
    mrc_fld_set_param_int_array(fld, "offs", 4,
			       (int[4]) { ilo[0], ilo[1], ilo[2], 0 });
    mrc_fld_set_param_int_array(fld, "dims", 4,
			       (int [4]) { ihi[0] - ilo[0], ihi[1] - ilo[1], ihi[2] - ilo[2],
				    mrc_fld_nr_comps(gfld) });
    mrc_fld_setup(fld);

    MPI_Datatype dtype;
    switch (mrc_fld_data_type(fld)) {
      case MRC_NT_FLOAT:
       dtype = MPI_FLOAT;
        for (int m = 0; m < mrc_fld_nr_comps(gfld); m++) {
          for (int iz = ilo[2]; iz < ihi[2]; iz++) {
            for (int iy = ilo[1]; iy < ihi[1]; iy++) {
              for (int ix = ilo[0]; ix < ihi[0]; ix++) {
                MRC_S4(fld,ix,iy,iz,m) = MRC_S4(gfld,ix,iy,iz,m);
              }
            }
          }
        }
      break;
      case MRC_NT_DOUBLE: 
        dtype = MPI_DOUBLE;
        for (int m = 0; m < mrc_fld_nr_comps(gfld); m++) {
          for (int iz = ilo[2]; iz < ihi[2]; iz++) {
            for (int iy = ilo[1]; iy < ihi[1]; iy++) {
              for (int ix = ilo[0]; ix < ihi[0]; ix++) {
                MRC_D4(fld,ix,iy,iz,m) = MRC_D4(gfld,ix,iy,iz,m);
              }
            }
          }
        }
        break;
      case MRC_NT_INT:
        dtype = MPI_INT;
        for (int m = 0; m < mrc_fld_nr_comps(gfld); m++) {
          for (int iz = ilo[2]; iz < ihi[2]; iz++) {
            for (int iy = ilo[1]; iy < ihi[1]; iy++) {
              for (int ix = ilo[0]; ix < ihi[0]; ix++) {
                MRC_I4(fld,ix,iy,iz,m) = MRC_I4(gfld,ix,iy,iz,m);
              }
            }
          }
        }
        break;
      default: assert(0);
    }
    
    MPI_Isend(fld->_nd->arr, mrc_fld_len(fld), dtype, block->rank, block->patch,
	      mrc_io_comm(io), &send->reqs[i]);
    block->fld = fld;
  }
}

static void
collective_m3_send_end(struct mrc_io *io, struct collective_m3_ctx *ctx)
{
  struct mrc_redist_read *send = &ctx->read_send;

  MPI_Waitall(send->n, send->reqs, MPI_STATUSES_IGNORE);
  for (int i = 0; i < send->n; i++) {
    mrc_fld_destroy(send->blocks[i].fld);
  }
  free(send->blocks);
  free(send->reqs);
}

// ----------------------------------------------------------------------

static void
collective_m3_recv_setup(struct mrc_io *io, struct collective_m3_ctx *ctx,
			 struct mrc_domain *domain)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct mrc_redist_read *recv = &ctx->read_recv;

  recv->n = 0;
  for (int p = 0; p < ctx->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(domain, p, &info);

    for (int writer = 0; writer < xdmf->nr_writers; writer++) {
      int writer_off[3], writer_dims[3];
      get_writer_off_dims(ctx, writer, writer_off, writer_dims);
      int ilo[3], ihi[3];
      if (find_intersection(ilo, ihi, info.off, info.ldims,
			    writer_off, writer_dims)) {
	recv->n++;
      }
    }
  }

  recv->reqs = calloc(recv->n, sizeof(*recv->reqs));
  recv->blocks = calloc(recv->n, sizeof(*recv->blocks));

  int i = 0;
  for (int p = 0; p < ctx->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(domain, p, &info);

    for (int writer = 0; writer < xdmf->nr_writers; writer++) {
      if (i == recv->n) {
	break;
      }
      struct collective_m3_entry *block = &recv->blocks[i];

      int writer_off[3], writer_dims[3];
      get_writer_off_dims(ctx, writer, writer_off, writer_dims);
      if (!find_intersection(block->ilo, block->ihi, info.off, info.ldims,
			     writer_off, writer_dims)) {
	continue;
      }
      block->rank = xdmf->writers[writer];
      block->patch = p;
      block->global_patch = info.global_patch;
      // patch-local indices from here on
      for (int d = 0; d < 3; d++) {
	block->ilo[d] -= info.off[d];
	block->ihi[d] -= info.off[d];
      }
      i++;
    }
  }
}

static void
collective_m3_recv_begin(struct mrc_io *io, struct collective_m3_ctx *ctx,
			 struct mrc_domain *domain, struct mrc_fld *m3)
{
  struct mrc_redist_read *recv = &ctx->read_recv;

  collective_m3_recv_setup(io, ctx, domain);

  for (int i = 0; i < recv->n; i++) {
    struct collective_m3_entry *block = &recv->blocks[i];

    struct mrc_fld *fld = mrc_fld_create(MPI_COMM_NULL);
    int *ilo = block->ilo, *ihi = block->ihi; // FIXME, -> off, dims
    mrc_fld_set_type(fld, mrc_fld_type(m3));
    mrc_fld_set_param_int_array(fld, "offs", 4,
			       (int[4]) { ilo[0], ilo[1], ilo[2], 0 });
    mrc_fld_set_param_int_array(fld, "dims", 4,
			       (int[4]) { ihi[0] - ilo[0], ihi[1] - ilo[1], ihi[2] - ilo[2],
				   mrc_fld_nr_comps(m3) });
    mrc_fld_setup(fld);

    MPI_Datatype dtype;
    switch (mrc_fld_data_type(fld)) {
      case MRC_NT_FLOAT: dtype = MPI_FLOAT; break;
      case MRC_NT_DOUBLE: dtype = MPI_DOUBLE; break;
      case MRC_NT_INT: dtype = MPI_INT; break;
      default: assert(0);
    }
    
    MPI_Irecv(fld->_nd->arr, mrc_fld_len(fld), dtype, block->rank,
	      block->global_patch, mrc_io_comm(io), &recv->reqs[i]);
    block->fld = fld;
  }
}

static void
collective_m3_recv_end(struct mrc_io *io, struct collective_m3_ctx *ctx,
		       struct mrc_domain *domain, struct mrc_fld *m3)
{
  struct mrc_redist_read *recv = &ctx->read_recv;

  MPI_Waitall(recv->n, recv->reqs, MPI_STATUSES_IGNORE);
  
  for (int i = 0; i < recv->n; i++) {
    struct collective_m3_entry *block = &recv->blocks[i];
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, block->patch);

    int *ilo = block->ilo, *ihi = block->ihi;
    switch (mrc_fld_data_type(m3)) {
      case MRC_NT_FLOAT:
        for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
          for (int iz = ilo[2]; iz < ihi[2]; iz++) {
            for (int iy = ilo[1]; iy < ihi[1]; iy++) {
              for (int ix = ilo[0]; ix < ihi[0]; ix++) {
                MRC_S5((m3p)->_fld, ix, iy, iz, m, (m3p)->_p) = MRC_S4(block->fld,ix,iy,iz,m);
                }
              }
            }
          }
        break;
      case MRC_NT_DOUBLE:
        for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
          for (int iz = ilo[2]; iz < ihi[2]; iz++) {
            for (int iy = ilo[1]; iy < ihi[1]; iy++) {
              for (int ix = ilo[0]; ix < ihi[0]; ix++) {
                MRC_D5((m3p)->_fld, ix, iy, iz, m, (m3p)->_p) = MRC_D4(block->fld,ix,iy,iz,m);
                }
              }
            }
          }
        break;
      case MRC_NT_INT:
        for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
          for (int iz = ilo[2]; iz < ihi[2]; iz++) {
            for (int iy = ilo[1]; iy < ihi[1]; iy++) {
              for (int ix = ilo[0]; ix < ihi[0]; ix++) {
                MRC_I5((m3p)->_fld, ix, iy, iz, m, (m3p)->_p) = MRC_I4(block->fld,ix,iy,iz,m);
                }
              }
            }
          }
        break;
      }
    mrc_fld_patch_put(m3);
    mrc_fld_destroy(block->fld);
  }
  free(recv->blocks);
  free(recv->reqs);
}

struct read_m3_cb_data {
  struct mrc_io *io;
  struct mrc_fld *gfld;
  hid_t filespace;
  hid_t memspace;
  hid_t dxpl;
};

static herr_t
read_m3_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{
  struct read_m3_cb_data *data = op_data;
  int ierr;

  hid_t group_fld = H5Gopen(g_id, name, H5P_DEFAULT); H5_CHK(group_fld);
  int m;
  ierr = H5LTget_attribute_int(group_fld, ".", "m", &m); CE;
  hid_t group = H5Gopen(group_fld, "p0", H5P_DEFAULT); H5_CHK(group);

  hid_t dset = H5Dopen(group, "3d", H5P_DEFAULT); H5_CHK(dset);
  struct mrc_fld *gfld = data->gfld;
  int *ib = gfld->_ghost_offs;
  hid_t dtype;
  void *arr;  
  switch (mrc_fld_data_type(gfld)) {
  case MRC_NT_FLOAT:
    dtype = H5T_NATIVE_FLOAT;
    arr = &MRC_S4(gfld, ib[0],ib[1],ib[2],m);
    break;
  case MRC_NT_DOUBLE:
    dtype = H5T_NATIVE_DOUBLE;
    arr = &MRC_D4(gfld, ib[0],ib[1],ib[2],m); 
    break;
  case MRC_NT_INT:
    dtype = H5T_NATIVE_INT;
    arr = &MRC_I4(gfld, ib[0],ib[1],ib[2],m);
    break;
  default: assert(0);
  }
  assert(!gfld->_aos);

  ierr = H5Dread(dset, dtype, data->memspace, data->filespace,
		 data->dxpl, arr); CE;
  ierr = H5Dclose(dset); CE;
  
  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group_fld); CE;

  return 0;
}

// ----------------------------------------------------------------------
// xdmf_collective_read_m3

static void
collective_m3_read_fld(struct mrc_io *io, struct collective_m3_ctx *ctx,
		       hid_t group0, struct mrc_fld *fld)
{
  struct xdmf *xdmf = to_xdmf(io);
  int ierr;

  int writer_rank;
  MPI_Comm_rank(xdmf->comm_writers, &writer_rank);
  int writer_dims[3], writer_off[3];
  get_writer_off_dims(ctx, writer_rank, writer_off, writer_dims);
  /* mprintf("writer_off %d %d %d dims %d %d %d\n", */
  /* 	    writer_off[0], writer_off[1], writer_off[2], */
  /* 	    writer_dims[0], writer_dims[1], writer_dims[2]); */

  int nr_comps = mrc_fld_nr_comps(fld);
  mrc_fld_set_param_int_array(fld, "dims", 4,
			      (int[4]) { writer_dims[0], writer_dims[1], writer_dims[2], nr_comps });
  mrc_fld_set_param_int_array(fld, "offs", 4,
			      (int[4]) { writer_off[0], writer_off[1], writer_off[2], 0 });
  mrc_fld_setup(fld);

  hsize_t fdims[3] = { ctx->gdims[2], ctx->gdims[1], ctx->gdims[0] };
  hsize_t foff[3] = { writer_off[2], writer_off[1], writer_off[0] };
  hsize_t mdims[3] = { writer_dims[2], writer_dims[1], writer_dims[0] };
  hid_t filespace = H5Screate_simple(3, fdims, NULL); H5_CHK(filespace);
  ierr = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL,
			     mdims, NULL); CE;
  hid_t memspace = H5Screate_simple(3, mdims, NULL); H5_CHK(memspace);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl);

#ifdef H5_HAVE_PARALLEL
  if (xdmf->use_independent_io) {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
  } else {
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  }
#endif

  struct read_m3_cb_data cb_data = {
    .io        = io,
    .gfld      = fld,
    .memspace  = memspace,
    .filespace = filespace,
    .dxpl      = dxpl,
  };
  
  hsize_t idx = 0;
  H5Literate_by_name(group0, ".", H5_INDEX_NAME, H5_ITER_INC, &idx,
		     read_m3_cb, &cb_data, H5P_DEFAULT);
  
  ierr = H5Pclose(dxpl); CE;
  ierr = H5Sclose(memspace); CE;
  ierr = H5Sclose(filespace); CE;
}

static void
xdmf_collective_read_m3(struct mrc_io *io, const char *path, struct mrc_fld *m3)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  //assert(m3->_data_type == MRC_NT_FLOAT);
  struct collective_m3_ctx ctx;
  collective_m3_init(io, &ctx, m3->_domain);

  // If we have an aos field, we need to get it as soa for reading and distribution
  struct mrc_fld *m3_soa = m3;
  if (m3->_aos) {
    switch (mrc_fld_data_type(m3)) {
      case MRC_NT_FLOAT: m3_soa = mrc_fld_get_as(m3, "float"); break;
      case MRC_NT_DOUBLE: m3_soa = mrc_fld_get_as(m3, "double"); break;
      case MRC_NT_INT: m3_soa = mrc_fld_get_as(m3, "int"); break;
      default: assert(0);
    }
  }

  if (xdmf->is_writer) {
    struct mrc_fld *gfld = mrc_fld_create(MPI_COMM_SELF);
    mrc_fld_set_param_int(gfld, "nr_comps", mrc_fld_nr_comps(m3));
    switch (mrc_fld_data_type(m3)) {
      case MRC_NT_FLOAT: mrc_fld_set_type(gfld, "float"); break;
      case MRC_NT_DOUBLE: mrc_fld_set_type(gfld, "double"); break;
      case MRC_NT_INT: mrc_fld_set_type(gfld, "int"); break;
      default: assert(0);
    }

    hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group0);
    collective_m3_read_fld(io, &ctx, group0, gfld);
    ierr = H5Gclose(group0); CE;

    collective_m3_recv_begin(io, &ctx, m3->_domain, m3_soa);
    collective_m3_send_begin(io, &ctx, m3->_domain, gfld);
    collective_m3_recv_end(io, &ctx, m3->_domain, m3_soa);
    collective_m3_send_end(io, &ctx);
    mrc_fld_destroy(gfld);
  } else {
    collective_m3_recv_begin(io, &ctx, m3->_domain, m3_soa);
    collective_m3_recv_end(io, &ctx, m3->_domain, m3_soa);
  }

  if (m3->_aos) {
    mrc_fld_put_as(m3_soa, m3);
  }
}

// ======================================================================
// mrc_io_ops_xdmf_collective

struct mrc_io_ops mrc_io_xdmf_collective_ops = {
  .name          = "xdmf_collective",
  .size          = sizeof(struct xdmf),
  .param_descr   = xdmf_collective_descr,
  .parallel      = true,
  .setup         = xdmf_collective_setup,
  .destroy       = xdmf_collective_destroy,
  .open          = xdmf_collective_open,
  .close         = xdmf_collective_close,
  .write_attr    = xdmf_collective_write_attr,
  .read_attr     = xdmf_collective_read_attr,
  .write_m1      = xdmf_collective_write_m1,
  .read_m1       = xdmf_collective_read_m1,
  .write_m3      = xdmf_collective_write_m3,
  .read_m3       = xdmf_collective_read_m3,
};


