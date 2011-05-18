
#include <mrc_io_private.h>
#include <mrc_params.h>
#include "mrc_io_xdmf_lib.h"

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <string.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

struct xdmf_file {
  hid_t h5_file;
  list_t xdmf_spatial_list;
};

struct xdmf {
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
  {},
};
#undef VAR

#define to_xdmf(io) mrc_to_subobj(io, struct xdmf)

// ----------------------------------------------------------------------
// xdmf_collective_setup

static void
xdmf_collective_setup(struct mrc_io *io)
{
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
    ierr = H5LTset_attribute_int(group, ".", name, &pv->u_int, 1); CE;
    break;
  case PT_BOOL: {
    int val = pv->u_bool;
    ierr = H5LTset_attribute_int(group, ".", name, &val, 1); CE;
    break;
  }
  case PT_FLOAT:
    ierr = H5LTset_attribute_float(group, ".", name, &pv->u_float, 1); CE;
    break;
  case PT_DOUBLE:
    ierr = H5LTset_attribute_double(group, ".", name, &pv->u_double, 1); CE;
    break;
  case PT_STRING:
    ierr = H5LTset_attribute_string(group, ".", name, pv->u_string); CE;
    break;
  case PT_INT3:
    ierr = H5LTset_attribute_int(group, ".", name, pv->u_int3, 3); CE;
    break;
  case PT_FLOAT3:
    ierr = H5LTset_attribute_float(group, ".", name, pv->u_float3, 3); CE;
    break;
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
      ierr = H5LTget_attribute_int(group, ".", name, &pv->u_int); CE;
      break;
    case PT_BOOL: ;
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
      break;
    case PT_INT3:
      ierr = H5LTget_attribute_int(group, ".", name, pv->u_int3); CE;
      break;
    case PT_FLOAT3:
      ierr = H5LTget_attribute_float(group, ".", name, pv->u_float3); CE;
      break;
    }
    ierr = H5Gclose(group); CE;
  }

  int root = xdmf->writers[0];
  MPI_Comm comm = mrc_io_comm(io);
  switch (type) {
  case PT_SELECT:
  case PT_INT:
    MPI_Bcast(&pv->u_int, 1, MPI_INT, root, comm);
    break;
  case PT_BOOL: ;
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
      len = strlen(pv->u_string);
    }
    MPI_Bcast(&len, 1, MPI_INT, root, comm);
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
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// collective_m1_write_f1
// does the actual write of the f1 to the file
// only called on writer procs

static void
collective_m1_write_f1(struct mrc_io *io, const char *path, struct mrc_f1 *f1,
		       int m, hid_t group0)
{
  int ierr;

  hid_t group_fld = H5Gcreate(group0, mrc_f1_comp_name(f1, m), H5P_DEFAULT,
			      H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);
  ierr = H5LTset_attribute_int(group_fld, ".", "m", &m, 1); CE;
  
  hid_t group = H5Gcreate(group_fld, "p0", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
  int i0 = 0;
  ierr = H5LTset_attribute_int(group, ".", "global_patch", &i0, 1); CE;

  hsize_t fdims[1] = { mrc_f1_ghost_dims(f1)[0] };
  hid_t filespace = H5Screate_simple(1, fdims, NULL); H5_CHK(filespace);
  hid_t dset = H5Dcreate(group, "1d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl); // FIXME, consolidate
#ifdef H5_HAVE_PARALLEL
  struct xdmf *xdmf = to_xdmf(io);
  if (xdmf->use_independent_io) {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT); CE;
  } else {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); CE;
  }
#endif
  ierr = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, dxpl, f1->arr); CE;
  
  ierr = H5Dclose(dset); CE;
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
};

static void
collective_m1_send_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			 struct mrc_m1 *m1, int m)
{
  struct xdmf *xdmf = to_xdmf(io);
  int dim = ctx->dim;

  ctx->send_reqs = calloc(ctx->nr_patches, sizeof(*ctx->send_reqs));
  ctx->nr_send_reqs = 0;

  for (int p = 0; p < ctx->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(m1->domain, p, &info);
    bool skip = false;
    for (int d = 0; d < 3; d++) {
      if (d != dim && info.off[d] != 0) {
	skip = true;
      }
    }
    if (skip) {
      continue;
    }
    
    struct mrc_m1_patch *m1p = mrc_m1_patch_get(m1, p);
    int ib = 0;
    if (info.off[dim] == 0) { // FIXME, -> generic code
      ib -= ctx->sw;
    }
    int ie = info.ldims[dim];
    if (info.off[dim] + info.ldims[dim] == ctx->gdims[dim]) {
      ie += ctx->sw;
    }
    //    mprintf("send to %d tag %d\n", xdmf->writers[0], info.global_patch);
    MPI_Isend(&MRC_M1(m1p, m, ib), ie - ib, MPI_FLOAT,
	      xdmf->writers[0], info.global_patch, mrc_io_comm(io),
	      &ctx->send_reqs[ctx->nr_send_reqs++]);
    mrc_m1_patch_put(m1);
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
			 struct mrc_domain *domain, struct mrc_f1 *f1, int m)
{
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
      ib -= ctx->sw;
    }
    int ie = info.off[dim] + info.ldims[dim];
    if (ie == ctx->gdims[dim]) {
      ie += ctx->sw;
    }
    //    mprintf("recv from %d tag %d\n", info.rank, gp);
    MPI_Irecv(&MRC_F1(f1, 0, ib), ie - ib, MPI_FLOAT, info.rank,
	      gp, mrc_io_comm(io), &ctx->recv_reqs[ctx->nr_recv_reqs++]);
  }
  assert(ctx->nr_recv_reqs == ctx->np[dim]);
}

static void
collective_m1_recv_end(struct mrc_io *io, struct collective_m1_ctx *ctx)
{
  MPI_Waitall(ctx->nr_recv_reqs, ctx->recv_reqs, MPI_STATUSES_IGNORE);
  free(ctx->recv_reqs);
}

// ----------------------------------------------------------------------
// xdmf_collective_write_m1

static void
xdmf_collective_write_m1(struct mrc_io *io, const char *path, struct mrc_m1 *m1)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  if (xdmf->nr_writers > 1) {
    MHERE; // FIXME
    return;
  }
  assert(xdmf->nr_writers == 1);
  struct collective_m1_ctx ctx;
  int nr_comps;
  mrc_m1_get_param_int(m1, "nr_comps", &nr_comps);
  mrc_m1_get_param_int(m1, "dim", &ctx.dim);
  mrc_m1_get_param_int(m1, "sw", &ctx.sw);
  mrc_domain_get_global_dims(m1->domain, ctx.gdims);
  mrc_domain_get_nr_global_patches(m1->domain, &ctx.nr_global_patches);
  mrc_domain_get_nr_procs(m1->domain, ctx.np);
  mrc_domain_get_patches(m1->domain, &ctx.nr_patches);

  int dim = ctx.dim;
  if (xdmf->is_writer) {
    struct mrc_f1 *f1 = mrc_f1_create(MPI_COMM_SELF);
    mrc_f1_set_param_int(f1, "dimsx", ctx.gdims[dim]);
    mrc_f1_set_param_int(f1, "sw", ctx.sw);
    mrc_f1_setup(f1);

    hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group0);
    for (int m = 0; m < nr_comps; m++) {
      mrc_f1_set_comp_name(f1, 0, mrc_m1_comp_name(m1, m));

      collective_m1_recv_begin(io, &ctx, m1->domain, f1, m);
      collective_m1_send_begin(io, &ctx, m1, m);
      collective_m1_recv_end(io, &ctx);
      collective_m1_write_f1(io, path, f1, m, group0);
      collective_m1_send_end(io, &ctx);
    }
    ierr = H5Gclose(group0); CE;

    mrc_f1_destroy(f1);
  } else { // not writer
    for (int m = 0; m < nr_comps; m++) {
      collective_m1_send_begin(io, &ctx, m1, m);
      collective_m1_send_end(io, &ctx);
    }
  }
}

// ----------------------------------------------------------------------

static void
collective_m1_read_recv_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			      struct mrc_m1 *m1, int m)
{
  struct xdmf *xdmf = to_xdmf(io);

  ctx->recv_reqs = calloc(ctx->nr_patches, sizeof(*ctx->recv_reqs));
  ctx->nr_recv_reqs = 0;
  for (int p = 0; p < ctx->nr_patches; p++) {
    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(m1->domain, p, &info);
    int ib = -ctx->sw;
    int ie = info.ldims[ctx->dim] + ctx->sw;
    struct mrc_m1_patch *m1p = mrc_m1_patch_get(m1, p);
    //	mprintf("recv to %d tag %d\n", xdmf->writers[0], info.global_patch);
    MPI_Irecv(&MRC_M1(m1p, m, ib), ie - ib, MPI_FLOAT,
	      xdmf->writers[0], info.global_patch, mrc_io_comm(io),
	      &ctx->recv_reqs[ctx->nr_recv_reqs++]);
    mrc_m1_patch_put(m1);
  }
}

static void
collective_m1_read_send_begin(struct mrc_io *io, struct collective_m1_ctx *ctx,
			      struct mrc_domain *domain, struct mrc_f1 *f1, int m)
{
  int dim = ctx->dim;
  ctx->send_reqs = calloc(ctx->nr_global_patches, sizeof(*ctx->send_reqs));
  ctx->nr_send_reqs = 0;
  for (int gp = 0; gp < ctx->nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(domain, gp, &info);
    int ib = info.off[dim] - ctx->sw;
    int ie = info.off[dim] + info.ldims[dim] + ctx->sw;
    //	mprintf("send from %d tag %d\n", info.rank, gp);
    MPI_Isend(&MRC_F1(f1, m, ib), ie - ib, MPI_FLOAT,
	      info.rank, gp, mrc_io_comm(io), &ctx->send_reqs[ctx->nr_send_reqs++]);
  }
}

struct read_m1_cb_data {
  struct mrc_io *io;
  struct mrc_f1 *gfld;
  hid_t filespace;
  hid_t memspace;
  hid_t dxpl;
};

static herr_t
read_m1_cb(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data)
{
  struct read_m1_cb_data *data = op_data;
  int ierr;

  hid_t group_fld = H5Gopen(g_id, name, H5P_DEFAULT); H5_CHK(group_fld);
  int m;
  ierr = H5LTget_attribute_int(group_fld, ".", "m", &m); CE;
  hid_t group = H5Gopen(group_fld, "p0", H5P_DEFAULT); H5_CHK(group);

  hid_t dset = H5Dopen(group, "1d", H5P_DEFAULT); H5_CHK(dset);
  struct mrc_f1 *gfld = data->gfld;
  int *ib = gfld->_ghost_off;
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
xdmf_collective_read_m1(struct mrc_io *io, const char *path, struct mrc_m1 *m1)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  if (xdmf->nr_writers > 1) {
    MHERE; // FIXME
    return;
  }
  assert(xdmf->nr_writers == 1); // FIXME
  struct collective_m1_ctx ctx;
  int nr_comps, gdims[3];
  mrc_m1_get_param_int(m1, "nr_comps", &nr_comps);
  mrc_m1_get_param_int(m1, "dim", &ctx.dim);
  mrc_m1_get_param_int(m1, "sw", &ctx.sw);
  mrc_domain_get_global_dims(m1->domain, gdims);
  mrc_domain_get_nr_global_patches(m1->domain, &ctx.nr_global_patches);
  mrc_domain_get_patches(m1->domain, &ctx.nr_patches);

  if (xdmf->is_writer) {
    struct mrc_f1 *f1 = mrc_f1_create(MPI_COMM_SELF);
    mrc_f1_set_param_int(f1, "nr_comps", nr_comps);
    mrc_f1_set_param_int(f1, "dimsx", gdims[ctx.dim]);
    mrc_f1_set_param_int(f1, "sw", ctx.sw);
    mrc_f1_setup(f1);

    hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group0);

    hsize_t hgdims[1] = { mrc_f1_ghost_dims(f1)[0] };

    hid_t filespace = H5Screate_simple(1, hgdims, NULL); H5_CHK(filespace);
    hid_t memspace = H5Screate_simple(1, hgdims, NULL); H5_CHK(memspace);
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
      collective_m1_read_send_begin(io, &ctx, m1->domain, f1, m);
      collective_m1_recv_end(io, &ctx);
      collective_m1_send_end(io, &ctx);
    }
  } else { // not writer
    for (int m = 0; m < nr_comps; m++) {
      collective_m1_read_recv_begin(io, &ctx, m1, m);
      collective_m1_recv_end(io, &ctx);
    }
  }
}

// ======================================================================

// ----------------------------------------------------------------------
// collective_write_f3
// does the actual write of the partial f3 to the file
// only called on writer procs

static void
collective_write_f3(struct mrc_io *io, const char *path, struct mrc_f3 *f3, int m,
		    struct mrc_m3 *m3, struct xdmf_spatial *xs, hid_t group0)
{
  int ierr;
  int gdims[3];
  mrc_domain_get_global_dims(m3->domain, gdims);

  xdmf_spatial_save_fld_info(xs, strdup(m3->name[m]), strdup(path), false);

  hid_t group_fld = H5Gcreate(group0, m3->name[m], H5P_DEFAULT,
			      H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group_fld);
  ierr = H5LTset_attribute_int(group_fld, ".", "m", &m, 1); CE;
  
  hid_t group = H5Gcreate(group_fld, "p0", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
  int i0 = 0;
  ierr = H5LTset_attribute_int(group, ".", "global_patch", &i0, 1); CE;

  hsize_t fdims[3] = { gdims[2], gdims[1], gdims[0] };
  hid_t filespace = H5Screate_simple(3, fdims, NULL); H5_CHK(filespace);
  hid_t dset = H5Dcreate(group, "3d", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
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
  const int *im = mrc_f3_ghost_dims(f3), *ib = mrc_f3_ghost_off(f3);
  hsize_t mdims[3] = { im[2], im[1], im[0] };
  hsize_t foff[3] = { ib[2], ib[1], ib[0] };
  hid_t memspace = H5Screate_simple(3, mdims, NULL);
  ierr = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL,
			     mdims, NULL); CE;

  ierr = H5Dwrite(dset, H5T_NATIVE_FLOAT, memspace, filespace, dxpl, f3->arr); CE;
  
  ierr = H5Dclose(dset); CE;
  ierr = H5Sclose(memspace); CE;
  ierr = H5Sclose(filespace); CE;
  ierr = H5Pclose(dxpl); CE;

  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group_fld); CE;
}

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

static bool
find_intersection(int *ilo, int *ihi, const int *ib1, const int *im1,
		  const int *ib2, const int *im2)
{
  for (int d = 0; d < 3; d++) {
    ilo[d] = MAX(ib1[d], ib2[d]);
    ihi[d] = MIN(ib1[d] + im1[d], ib2[d] + im2[d]);
    if (ihi[d] - ilo[d] <= 0) {
      return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------
// collective helper context

struct collective_m3_ctx {
  int gdims[3];
  int writer_dims[3], writer_off[3];
  int nr_patches, nr_global_patches;
  int slow_dim;
  int slow_indices_per_writer;
  int slow_indices_rmndr;
  int *recv_gps;
  struct mrc_f3 **recv_f3s;
  int total_sends;
  MPI_Request *send_reqs;
  int total_recvs;
  MPI_Request *recv_reqs;
};

static void
get_writer_off_dims(struct collective_m3_ctx *ctx, int writer,
		    int *writer_off, int *writer_dims)
{
  for (int d = 0; d < 3; d++) {
    writer_dims[d] = ctx->gdims[d];
    writer_off[d] = 0;
  }
  writer_dims[ctx->slow_dim] = ctx->slow_indices_per_writer + (writer < ctx->slow_indices_rmndr);
  if (writer < ctx->slow_indices_rmndr) {
    writer_off[ctx->slow_dim] = (ctx->slow_indices_per_writer + 1) * writer;
  } else {
    writer_off[ctx->slow_dim] = ctx->slow_indices_rmndr +
      ctx->slow_indices_per_writer * writer;
  }
}

// ----------------------------------------------------------------------
// collective_send_f3_begin

static void
collective_send_f3_begin(struct collective_m3_ctx *ctx, struct mrc_io *io,
			 struct mrc_m3 *m3, int m)
{
  struct xdmf *xdmf = to_xdmf(io);

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m3->domain, &nr_patches);
  ctx->total_sends = 0;
  for (int p = 0; p < nr_patches; p++) {
    int *off = patches[p].off, *ldims = patches[p].ldims;
    for (int writer = 0; writer < xdmf->nr_writers; writer++) {
      // don't send to self
      if (xdmf->writers[writer] == io->rank) {
	continue;
      }
      int writer_off[3], writer_dims[3];
      get_writer_off_dims(ctx, writer, writer_off, writer_dims);
      int ilo[3], ihi[3];
      bool has_intersection =
	find_intersection(ilo, ihi, off, ldims, writer_off, writer_dims);
      if (has_intersection) {
	ctx->total_sends++;
      }
    }
  }
  //  mprintf("total_sends = %d\n", total_sends);

  int sr = 0;
  ctx->send_reqs = calloc(ctx->total_sends, sizeof(*ctx->send_reqs));

  for (int p = 0; p < nr_patches; p++) {
    int *off = patches[p].off, *ldims = patches[p].ldims;
    for (int writer = 0; writer < xdmf->nr_writers; writer++) {
      // don't send to self
      if (xdmf->writers[writer] == io->rank) {
	continue;
      }
      int writer_off[3], writer_dims[3];
      get_writer_off_dims(ctx, writer, writer_off, writer_dims);
      int ilo[3], ihi[3];
      bool has_intersection =
	find_intersection(ilo, ihi, off, ldims, writer_off, writer_dims);
      if (!has_intersection)
	continue;

      struct mrc_patch_info info;
      mrc_domain_get_local_patch_info(m3->domain, p, &info);
      struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
      /* mprintf("MPI_Isend -> %d gp %d len %d\n", xdmf->writers[writer], info.global_patch, */
      /* 	      m3p->im[0] * m3p->im[1] * m3p->im[2]); */
      MPI_Isend(&MRC_M3(m3p, m, m3p->ib[0], m3p->ib[1], m3p->ib[2]),
		m3p->im[0] * m3p->im[1] * m3p->im[2], MPI_FLOAT,
		xdmf->writers[writer], info.global_patch,
		mrc_io_comm(io), &ctx->send_reqs[sr++]);
      mrc_m3_patch_put(m3);
    }
  }
}

// ----------------------------------------------------------------------
// collective_send_f3_end

static void
collective_send_f3_end(struct collective_m3_ctx *ctx, struct mrc_io *io,
		       struct mrc_m3 *m3, int m)
{
  MPI_Waitall(ctx->total_sends, ctx->send_reqs, MPI_STATUSES_IGNORE);
  free(ctx->send_reqs);
}
    
// ----------------------------------------------------------------------
// collective_recv_f3_begin

static void
collective_recv_f3_begin(struct collective_m3_ctx *ctx,
			 struct mrc_io *io, struct mrc_f3 *f3,
			 struct mrc_m3 *m3)
{
  // find out who's sending, OPT: this way is not very scalable
  // could also be optimized by just looking at slow_dim
  // FIXME, figure out pattern and cache, at least across components

  int nr_global_patches;
  mrc_domain_get_nr_global_patches(m3->domain, &nr_global_patches);
  ctx->total_recvs = 0;
  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(m3->domain, gp, &info);
    // skip local patches for now
    if (info.rank == io->rank) {
      continue;
    }
    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     mrc_f3_ghost_off(f3), mrc_f3_ghost_dims(f3));
    if (has_intersection) {
      ctx->total_recvs++;
    }
  }

  //  mprintf("total_recvs = %d\n", total_recvs);
  int rr = 0;
  ctx->recv_gps = calloc(ctx->total_recvs, sizeof(*ctx->recv_gps));
  ctx->recv_reqs = calloc(ctx->total_recvs, sizeof(*ctx->recv_reqs));
  ctx->recv_f3s = calloc(ctx->total_recvs, sizeof(*ctx->recv_f3s));

  for (int gp = 0; gp < nr_global_patches; gp++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(m3->domain, gp, &info);
    // skip local patches for now
    if (info.rank == io->rank) {
      continue;
    }
    int ilo[3], ihi[3];
    int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					     mrc_f3_ghost_off(f3), mrc_f3_ghost_dims(f3));
    if (!has_intersection) {
      continue;
    }
    ctx->recv_gps[rr] = gp;
    struct mrc_f3 *recv_f3 = mrc_f3_create(MPI_COMM_NULL);
    mrc_f3_set_param_int3(recv_f3, "dims", info.ldims);
    mrc_f3_set_param_int(recv_f3, "sw", m3->sw);
    mrc_f3_setup(recv_f3);
    ctx->recv_f3s[rr] = recv_f3;
    
    /* mprintf("MPI_Irecv <- %d gp %d len %d\n", info.rank, info.global_patch, recv_f3->len); */
    MPI_Irecv(recv_f3->arr, recv_f3->len, MPI_FLOAT, info.rank,
	      info.global_patch, mrc_io_comm(io), &ctx->recv_reqs[rr++]);
  }
}

// ----------------------------------------------------------------------
// collective_recv_f3

static void
collective_recv_f3_end(struct collective_m3_ctx *ctx,
		       struct mrc_io *io, struct mrc_f3 *f3,
		       struct mrc_m3 *m3, int m)
{
  MPI_Waitall(ctx->total_recvs, ctx->recv_reqs, MPI_STATUSES_IGNORE);

  for (int rr = 0; rr < ctx->total_recvs; rr++) {
    struct mrc_patch_info info;
    mrc_domain_get_global_patch_info(m3->domain, ctx->recv_gps[rr], &info);
    int *off = info.off;
    // OPT, could be cached 2nd(?) and 3rd time
    int ilo[3], ihi[3];
    find_intersection(ilo, ihi, info.off, info.ldims,
		      mrc_f3_ghost_off(f3), mrc_f3_ghost_dims(f3));

    struct mrc_f3 *recv_f3 = ctx->recv_f3s[rr];
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_F3(f3,0, ix,iy,iz) =
	    MRC_F3(recv_f3, 0, ix - off[0], iy - off[1], iz - off[2]);
	}
      }
    }
    mrc_f3_destroy(recv_f3);
  }

  free(ctx->recv_reqs);
  free(ctx->recv_f3s);
  free(ctx->recv_gps);
}

// ----------------------------------------------------------------------
// collective_recv_f3_local

static void
collective_recv_f3_local(struct collective_m3_ctx *ctx,
			 struct mrc_io *io, struct mrc_f3 *f3,
			 struct mrc_m3 *m3, int m)
{
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m3->domain, &nr_patches);

  for (int p = 0; p < nr_patches; p++) {
    struct mrc_patch *patch = &patches[p];
    int *off = patch->off, *ldims = patch->ldims;

    int ilo[3], ihi[3];
    bool has_intersection =
      find_intersection(ilo, ihi, off, ldims, mrc_f3_ghost_off(f3), mrc_f3_ghost_dims(f3));
    if (!has_intersection) {
      continue;
    }
    struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
    for (int iz = ilo[2]; iz < ihi[2]; iz++) {
      for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	  MRC_F3(f3,0, ix,iy,iz) =
	    MRC_M3(m3p, m, ix - off[0], iy - off[1], iz - off[2]);
	}
      }
    }
    mrc_m3_patch_put(m3);
  }
}

// ----------------------------------------------------------------------
// xdmf_collective_write_m3

static void
xdmf_collective_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3)
{
  struct xdmf *xdmf = to_xdmf(io);

  struct collective_m3_ctx ctx;
  mrc_domain_get_global_dims(m3->domain, ctx.gdims);
  ctx.slow_dim = 2;
  while (ctx.gdims[ctx.slow_dim] == 1) {
    ctx.slow_dim--;
  }
  assert(ctx.slow_dim >= 0);
  int total_slow_indices = ctx.gdims[ctx.slow_dim];
  ctx.slow_indices_per_writer = total_slow_indices / xdmf->nr_writers;
  ctx.slow_indices_rmndr = total_slow_indices % xdmf->nr_writers;

  struct xdmf_file *file = &xdmf->file;
  struct xdmf_spatial *xs = xdmf_spatial_find(&file->xdmf_spatial_list,
					      mrc_domain_name(m3->domain));
  if (!xs) {
    xs = xdmf_spatial_create_m3_parallel(&file->xdmf_spatial_list,
					 mrc_domain_name(m3->domain),
					 m3->domain);
  }

  if (xdmf->is_writer) {
    int writer_rank;
    MPI_Comm_rank(xdmf->comm_writers, &writer_rank);
    int writer_dims[3], writer_off[3];
    get_writer_off_dims(&ctx, writer_rank, writer_off, writer_dims);
    /* mprintf("writer_off %d %d %d dims %d %d %d\n", */
    /* 	    writer_off[0], writer_off[1], writer_off[2], */
    /* 	    writer_dims[0], writer_dims[1], writer_dims[2]); */

    struct mrc_f3 *f3 = mrc_f3_create(MPI_COMM_NULL);
    mrc_f3_set_param_int3(f3, "off", writer_off);
    mrc_f3_set_param_int3(f3, "dims", writer_dims);
    mrc_f3_setup(f3);

    hid_t group0;
    if (H5Lexists(file->h5_file, path, H5P_DEFAULT) > 0) {
      group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT);
    } else {
      assert(0); // FIXME, can this happen?
      group0 = H5Gcreate(file->h5_file, path, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group0);
    }
    int nr_1 = 1;
    H5LTset_attribute_int(group0, ".", "nr_patches", &nr_1, 1);
    
    for (int m = 0; m < m3->nr_comp; m++) {
      collective_recv_f3_begin(&ctx, io, f3, m3);
      collective_send_f3_begin(&ctx, io, m3, m);
      collective_recv_f3_local(&ctx, io, f3, m3, m);
      collective_recv_f3_end(&ctx, io, f3, m3, m);
      collective_write_f3(io, path, f3, m, m3, xs, group0);
      collective_send_f3_end(&ctx, io, m3, m);
    }

    H5Gclose(group0);
    mrc_f3_destroy(f3);
  } else {
    for (int m = 0; m < m3->nr_comp; m++) {
      collective_send_f3_begin(&ctx, io, m3, m);
      collective_send_f3_end(&ctx, io, m3, m);
    }
  }
}

// ----------------------------------------------------------------------
// xdmf_collective_read_m3

struct read_m3_cb_data {
  struct mrc_io *io;
  struct mrc_f3 *gfld;
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
  struct mrc_f3 *gfld = data->gfld;
  int *ib = gfld->_ghost_off;
  ierr = H5Dread(dset, H5T_NATIVE_FLOAT, data->memspace, data->filespace,
		 data->dxpl, &MRC_F3(gfld, m, ib[0], ib[1], ib[2])); CE;
  ierr = H5Dclose(dset); CE;
  
  ierr = H5Gclose(group); CE;
  ierr = H5Gclose(group_fld); CE;

  return 0;
}

static void
collective_m3_read_f3(struct mrc_io *io, struct collective_m3_ctx *ctx,
		      hid_t group0, struct mrc_f3 *f3)
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

  mrc_f3_set_param_int3(f3, "off", writer_off);
  mrc_f3_set_param_int3(f3, "dims", writer_dims);
  mrc_f3_setup(f3);

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
    .gfld      = f3,
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
xdmf_collective_read_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3)
{
  struct xdmf *xdmf = to_xdmf(io);
  struct xdmf_file *file = &xdmf->file;
  int ierr;

  struct collective_m3_ctx ctx;
  // FIXME dupl
  mrc_domain_get_global_dims(m3->domain, ctx.gdims);
  mrc_domain_get_patches(m3->domain, &ctx.nr_patches);
  mrc_domain_get_nr_global_patches(m3->domain, &ctx.nr_global_patches);
  ctx.slow_dim = 2;
  while (ctx.gdims[ctx.slow_dim] == 1) {
    ctx.slow_dim--;
  }
  assert(ctx.slow_dim >= 0);
  int total_slow_indices = ctx.gdims[ctx.slow_dim];
  ctx.slow_indices_per_writer = total_slow_indices / xdmf->nr_writers;
  ctx.slow_indices_rmndr = total_slow_indices % xdmf->nr_writers;

  if (xdmf->is_writer) {
    struct mrc_f3 *gfld = mrc_f3_create(MPI_COMM_SELF);
    mrc_f3_set_param_int(gfld, "nr_comps", m3->nr_comp);

    hid_t group0 = H5Gopen(file->h5_file, path, H5P_DEFAULT); H5_CHK(group0);
    collective_m3_read_f3(io, &ctx, group0, gfld);
    int writer_rank;
    MPI_Comm_rank(xdmf->comm_writers, &writer_rank);
    get_writer_off_dims(&ctx, writer_rank, ctx.writer_off, ctx.writer_dims);

    ierr = H5Gclose(group0); CE;

    for (int gp = 0; gp < ctx.nr_global_patches; gp++) {
      struct mrc_patch_info info;
      mrc_domain_get_global_patch_info(m3->domain, gp, &info);

      int ilo[3], ihi[3];
      int has_intersection = find_intersection(ilo, ihi, info.off, info.ldims,
					       mrc_f3_ghost_off(gfld), mrc_f3_ghost_dims(gfld));
      if (!has_intersection) {
	continue;
      }
      mprintf("is %d:%d %d:%d %d:%d\n", ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);
      struct mrc_f3 *f3 = mrc_f3_create(MPI_COMM_NULL);
      mrc_f3_set_param_int3(f3, "off", ilo);
      mrc_f3_set_param_int3(f3, "dims",
			    (int [3]) { ihi[0] - ilo[0], ihi[1] - ilo[1], ihi[2] - ilo[2] });
      mrc_f3_set_param_int(f3, "nr_comps", m3->nr_comp);
      mrc_f3_setup(f3);
      assert(info.rank == io->rank);

      for (int m = 0; m < m3->nr_comp; m++) {
	for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	  for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	    for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	      MRC_F3(f3, m, ix,iy,iz) = MRC_F3(gfld, m, ix,iy,iz);
	    }
	  }
	}
      }

      struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, info.patch);
      int *off = info.off;
      for (int m = 0; m < m3->nr_comp; m++) {
	for (int iz = ilo[2]; iz < ihi[2]; iz++) {
	  for (int iy = ilo[1]; iy < ihi[1]; iy++) {
	    for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	      MRC_M3(m3p, m, ix-off[0],iy-off[1],iz-off[2]) = MRC_F3(f3, m, ix,iy,iz);
	    }
	  }
	}
      }
      mrc_m3_patch_put(m3);
    }

    mrc_f3_destroy(gfld);
  }
}

// ======================================================================
// mrc_io_ops_xdmf_collective

struct mrc_io_ops mrc_io_xdmf2_collective_ops = {
  .name          = "xdmf2_collective",
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


