
#include "mrc_io_private.h"
#include <mrc_params.h>
#include <mrc_profile.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

// TODO:
// 2D serial output (?)
// scaling for serial output
// vector/scalar

enum {
  ID_DIAGS = 14000,
  ID_DIAGS_CMD,
  ID_DIAGS_CREATE_OUTDIR,
  ID_DIAGS_CREATE_BASENAME,
  ID_DIAGS_CMD_CREATE,
  ID_DIAGS_CMD_DOMAIN_INFO,
  ID_DIAGS_CMD_OPEN,
  ID_DIAGS_CMD_WRITE,
  ID_DIAGS_CMD_WRITE_ATTR,
  ID_DIAGS_CMD_CRDX,
  ID_DIAGS_CMD_CRDY,
  ID_DIAGS_CMD_CRDZ,
  ID_DIAGS_BASENAME,
  ID_DIAGS_TIME,
  ID_DIAGS_FLDNAME,
  ID_DIAGS_SUBDOMAIN,
  ID_DIAGS_DATA,
  ID_DIAGS_2DDATA,
  ID_DIAGS_CMD_RESPONSE,
};

enum {
  DIAG_CMD_CREATE,
  DIAG_CMD_OPENFILE,
  DIAG_CMD_SHUTDOWN,
  DIAG_RESPONSE_SHUTDOWN_COMPLETE,
};

// ======================================================================
// diag client interface

struct diagc_combined_params {
  int rank_diagsrv;
};

#define VAR(x) (void *)offsetof(struct diagc_combined_params, x)

static struct param diagc_combined_params_descr[] = {
  { "rank_diagsrv"        , VAR(rank_diagsrv)      , PARAM_INT(0)       },
  {},
};

#undef VAR

// ----------------------------------------------------------------------
// diagc_combined_setup
//
// called once before first open, will do handshake with server

static void
diagc_combined_send_domain_info(struct mrc_io *io, struct mrc_domain *domain)
{
  if (io->diagc_domain_info_sent)
    return;

  struct diagc_combined_params *par = io->obj.subctx;

  int iw[9], *off = iw, *ldims = iw + 3, *gdims = iw + 6;
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(domain, &nr_patches);
  assert(nr_patches == 1);
  for (int d = 0; d < 3; d++) {
    off[d] = patches[0].off[d];
    ldims[d] = patches[0].ldims[d];
  }
  mrc_domain_get_global_dims(domain, gdims);
  if (io->rank == 0) {
    MPI_Send(iw, 0, MPI_CHAR, par->rank_diagsrv, ID_DIAGS_CMD_CREATE, MPI_COMM_WORLD);
  }
  MPI_Send(iw, 9, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD_DOMAIN_INFO, MPI_COMM_WORLD);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  for (int d = 0; d < 3; d++) {
    MPI_Send(&MRC_CRD(crds, d, 0), ldims[d], MPI_FLOAT, par->rank_diagsrv,
	     ID_DIAGS_CMD_CRDX + d, MPI_COMM_WORLD);
  }

  io->diagc_domain_info_sent = true;
}

static void
diagc_combined_setup(struct mrc_io *io)
{
  struct mrc_io_params *par_io = &io->par;
  struct diagc_combined_params *par = io->obj.subctx;

  mrc_io_setup_super(io);
  if (io->rank == 0) {
    int icmd[1] = { DIAG_CMD_CREATE };
    MPI_Send(icmd, 1, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD, MPI_COMM_WORLD);
    MPI_Send(par_io->outdir, strlen(par_io->outdir) + 1, MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_CREATE_OUTDIR, MPI_COMM_WORLD);
    MPI_Send(par_io->basename, strlen(par_io->basename) + 1, MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_CREATE_BASENAME, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_open
//
// opens a new output file (index itime)
// afterwards, diag_write_field() may be called repeatedly and
// must be concluded with a call to diag_close()

static void
diagc_combined_open(struct mrc_io *io, const char *mode)
{
  struct diagc_combined_params *par = io->obj.subctx;
  assert(strcmp(mode, "w") == 0); // only writing supported for now

  if (io->rank == 0) {
    int icmd[2] = { DIAG_CMD_OPENFILE, io->step };

    MPI_Send(icmd, 2, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD_OPEN, MPI_COMM_WORLD);
    MPI_Send(io->par.basename, strlen(io->par.basename) + 1, MPI_CHAR,
	     par->rank_diagsrv, ID_DIAGS_BASENAME, MPI_COMM_WORLD);
    MPI_Send(&io->time, 1, MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_TIME, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_close

static void
diagc_combined_close(struct mrc_io *io)
{
  struct diagc_combined_params *par = io->obj.subctx;

  if (io->rank == 0) {
    char str[] = "";
    MPI_Send(str, 1, MPI_CHAR, par->rank_diagsrv, ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_destroy
//
// shuts down the diag server process

static void
diagc_combined_destroy(struct mrc_io *io)
{
  static int pr_diag_shutdown = 0;
  int comm_world_rank = 0;
  int wait_for_response, icmd[2];

  struct diagc_combined_params *par = io->obj.subctx;

  MPI_Comm_rank(MPI_COMM_WORLD, &comm_world_rank);
  if (mrc_io_is_setup(io) && comm_world_rank == 0) {
    wait_for_response = 1;
  } else {
    wait_for_response = 0;
  }
  
  icmd[0] = DIAG_CMD_SHUTDOWN;
  icmd[1] = wait_for_response;
  MPI_Send(icmd, 2, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD, MPI_COMM_WORLD);

  if (0&&wait_for_response) {
    int response = 0;
    if (pr_diag_shutdown == 0) {
      pr_diag_shutdown = prof_register("diagc_combined_shutdown", 0, 0, 0);
    }
    prof_start(pr_diag_shutdown);
    MPI_Recv(&response, 1, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD_RESPONSE,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    assert(response == DIAG_RESPONSE_SHUTDOWN_COMPLETE);
    prof_stop(pr_diag_shutdown);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_write_field

static void
copy_and_scale(float *buf, struct mrc_fld *fld, int m, float scale)
{
  struct mrc_fld *f = mrc_fld_get_as(fld, "float");
  int i = 0;
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    buf[i++] = scale * MRC_F3(f, m, ix,iy,iz);
  } mrc_fld_foreach_end;
  mrc_fld_put_as(f, fld);
}

static void
diagc_combined_write_field(struct mrc_io *io, const char *path,
			   float scale, struct mrc_fld *fld, int m)
{
  struct diagc_combined_params *par = io->obj.subctx;

  diagc_combined_send_domain_info(io, fld->_domain);

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(fld->_domain, &nr_patches);
  assert(nr_patches == 1);
  int *ldims = patches[0].ldims;
  int nout = ldims[0] * ldims[1] * ldims[2];
  float *buf = calloc(sizeof(float), nout);

  if (io->rank == 0) {
    MPI_Send((char *)mrc_fld_name(fld), strlen(mrc_fld_name(fld)) + 1,
	     MPI_CHAR, par->rank_diagsrv, ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
    MPI_Send((char *)mrc_fld_comp_name(fld, m), strlen(mrc_fld_comp_name(fld, m)) + 1,
	     MPI_CHAR, par->rank_diagsrv, ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
    int outtype = DIAG_TYPE_3D;
    MPI_Send(&outtype, 1, MPI_INT, par->rank_diagsrv,
	     ID_DIAGS_CMD_WRITE, MPI_COMM_WORLD);
  }

  copy_and_scale(buf, fld, m, scale);

  int iw[6], *off = iw, *dims = iw + 3; // off, then dims
  for (int d = 0; d < 3; d++) {
    off[d] = patches[0].off[d];
    dims[d] = patches[0].ldims[d];
  }

  MPI_Send(iw, 6, MPI_INT, par->rank_diagsrv, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD);
  MPI_Send(buf, nout, MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_DATA, MPI_COMM_WORLD);

  free(buf);
}

static void
diagc_combined_write_m3(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  int nr_comps = mrc_fld_nr_comps(fld);
  for (int m = 0; m < nr_comps; m++) {
    assert(mrc_fld_comp_name(fld, m));
    diagc_combined_write_field(io, path, 1.f, fld, m);
  }
}

static void
diagc_combined_write_field2d(struct mrc_io *io, float scale, struct mrc_fld *fld,
			     int outtype, float sheet)
{
  struct diagc_combined_params *par = io->obj.subctx;

  diagc_combined_send_domain_info(io, fld->_domain);

  assert(outtype >= DIAG_TYPE_2D_X && outtype <= DIAG_TYPE_2D_Z);
  int dim = outtype - DIAG_TYPE_2D_X;

  if (io->rank == 0) {
    MPI_Send("mrc_f2", strlen("mrc_f2") + 1, MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
    MPI_Send((char *) mrc_fld_comp_name(fld, 0), strlen(mrc_fld_comp_name(fld, 0)) + 1,
	     MPI_CHAR, par->rank_diagsrv, ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
    MPI_Send(&outtype, 1, MPI_INT, par->rank_diagsrv,
	     ID_DIAGS_CMD_WRITE, MPI_COMM_WORLD);
    MPI_Send(&sheet, 1, MPI_FLOAT, par->rank_diagsrv,
	     ID_DIAGS_CMD_WRITE, MPI_COMM_WORLD);
  }

  int iw[6] = { -1, };
  if (mrc_fld_len(fld)) { // part of the slice?
    int *off = iw, *dims = iw + 3; // off, then dims
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(fld->_domain, &nr_patches);
    assert(nr_patches == 1);
    for (int d = 0; d < 3; d++) {
      off[d] = patches[0].off[d];
      dims[d] = patches[0].ldims[d];
    }
    off[dim] = 0;
    dims[dim] = 1;

    MPI_Send(iw, 6, MPI_INT, par->rank_diagsrv, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD);
    struct mrc_fld *f = mrc_fld_get_as(fld, "float");
    MPI_Send(fld->_nd->arr, mrc_fld_len(fld), MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_2DDATA, MPI_COMM_WORLD);
    mrc_fld_put_as(f, fld);
  } else {
    MPI_Send(iw, 6, MPI_INT, par->rank_diagsrv, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD);
  }
}

static void
diagc_combined_write_attr(struct mrc_io *io, const char *path, int type,
			  const char *name, union param_u *pv)
{
  struct diagc_combined_params *par = io->obj.subctx;

  if (io->rank == 0) {
    MPI_Send((char *)path, strlen(path) + 1, MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
    MPI_Send(&type, 1, MPI_INT, par->rank_diagsrv,
	     ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
    MPI_Send((char *)name, strlen(name) + 1, MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
    switch (type) {
    case PT_BOOL:
    case PT_INT:
    case PT_SELECT:
    case MRC_VAR_INT:
    case MRC_VAR_BOOL:
      MPI_Send(&pv->u_int, 1, MPI_INT, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_FLOAT:
      MPI_Send(&pv->u_float, 1, MPI_FLOAT, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_DOUBLE:
    case MRC_VAR_DOUBLE:
      MPI_Send(&pv->u_double, 1, MPI_DOUBLE, par->rank_diagsrv,
         ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_STRING:
      MPI_Send((char *)pv->u_string, strlen(pv->u_string) + 1, MPI_CHAR, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_INT3:
      MPI_Send(pv->u_int3, 3, MPI_INT, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_FLOAT3:
      MPI_Send(pv->u_float3, 3, MPI_FLOAT, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_DOUBLE3:
    case MRC_VAR_DOUBLE3:
      MPI_Send(pv->u_double3, 3, MPI_DOUBLE, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    case PT_INT_ARRAY:
      MPI_Send(&pv->u_int_array.nr_vals, 1, MPI_INT, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      MPI_Send(pv->u_int_array.vals, pv->u_int_array.nr_vals, MPI_INT, par->rank_diagsrv,
	       ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD);
      break;
    default:
      mprintf("type %d\n", type);
      assert(0);
    }
  }
}

struct mrc_io_ops mrc_io_combined_ops = {
  .name          = "combined",
  .size          = sizeof(struct diagc_combined_params),
  .param_descr   = diagc_combined_params_descr,
  .setup         = diagc_combined_setup,
  .destroy       = diagc_combined_destroy,
  .open          = diagc_combined_open,
  .write_m3      = diagc_combined_write_m3,
  .write_field2d = diagc_combined_write_field2d,
  .write_attr    = diagc_combined_write_attr,
  .close         = diagc_combined_close,
};

// ======================================================================
// diag server one

// generic diagnostics server running on one node.
// a particular implementation needs to take care of opening files/
// writing/closing them,
// the generic part here (in conjunction with the client side above)
// takes care of communicating / assembling the fields

// ----------------------------------------------------------------------

struct diagsrv_one {
  list_t mrc_io_list;

  // only valid from open() -> close()
  struct mrc_io *io;

  struct diagsrv_srv_ops *srv_ops;
  void *srv;
};

struct diagsrv_srv_ops {
  const char *name;
  void  (*create)(struct diagsrv_one *ds);
  void  (*set_domain)(struct diagsrv_one *ds, struct mrc_domain *domain);
  void  (*destroy)(struct diagsrv_one *ds);
  void  (*open)(struct diagsrv_one *ds, int step, float time);
  struct mrc_fld *(*get_gfld_2d)(struct diagsrv_one *ds, int dims[2]);
  struct mrc_fld *(*get_gfld_3d)(struct diagsrv_one *ds, int dims[3]);
  void  (*put_gfld_2d)(struct diagsrv_one *ds, struct mrc_fld *fld, char *fld_name,
		       int outtype, float sheet);
  void  (*put_gfld_3d)(struct diagsrv_one *ds, struct mrc_fld *fld);
  void  (*write_attr)(struct diagsrv_one *ds, const char *path, int type,
		      const char *name, union param_u *pv);
  void  (*close)(struct diagsrv_one *ds);
};

// ----------------------------------------------------------------------

struct mrc_io_entry {
  struct mrc_io *io;
  struct mrc_domain *domain;
  int ldims[3];

  char *basename;
  list_t entry;
};

static struct mrc_io_entry *
find_diagsrv_io(struct diagsrv_one *ds, const char *basename)
{
  struct mrc_io_entry *p;
  __list_for_each_entry(p, &ds->mrc_io_list, entry, struct mrc_io_entry) {
    if (strcmp(p->basename, basename) == 0) {
      return p;
    }
  }
  return NULL;
}

static void
create_diagsrv_io(struct diagsrv_one *ds,
		  const char *format, const char *outdir, const char *basename)
{
  struct mrc_io_entry *p = calloc(1, sizeof(*p));
  p->io = mrc_io_create(MPI_COMM_SELF);
  mrc_io_set_type(p->io, format);
  mrc_io_set_param_string(p->io, "outdir", outdir);
  mrc_io_set_param_string(p->io, "basename", basename);
  mrc_io_setup(p->io);
  mrc_io_view(p->io);

  p->basename = strdup(basename);
  list_add_tail(&p->entry, &ds->mrc_io_list);
}

// ----------------------------------------------------------------------
// diagsrv_one

struct diagsrv_srv {
  float *gfld;
  struct mrc_domain *domain;
};

static void
ds_srv_create(struct diagsrv_one *ds)
{
  struct diagsrv_srv *srv = malloc(sizeof(*srv));
  ds->srv = srv;
}

static void
ds_srv_set_domain(struct diagsrv_one *ds, struct mrc_domain *domain)
{
  struct diagsrv_srv *srv = ds->srv;
  srv->domain = domain;
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int nglobal = gdims[0] * gdims[1] * gdims[2];
  srv->gfld = malloc(nglobal * sizeof(float));
}

static void
ds_srv_destroy(struct diagsrv_one *ds)
{
  struct diagsrv_srv *srv = (struct diagsrv_srv *) ds->srv;
  free(srv->gfld);
  free(srv);
}

static void
ds_srv_open(struct diagsrv_one *ds, int step, float time)
{
  mrc_io_open(ds->io, "w", step, time);
}

static struct mrc_fld *
ds_srv_get_gfld_3d(struct diagsrv_one *ds, int gdims[3])
{
  struct diagsrv_srv *srv = (struct diagsrv_srv *) ds->srv;
  struct mrc_fld *fld = mrc_domain_fld_create(srv->domain, SW_0, NULL);
  mrc_fld_set_array(fld, srv->gfld);
  mrc_fld_setup(fld);
  return fld;
}

static struct mrc_fld *
ds_srv_get_gfld_2d(struct diagsrv_one *ds, int gdims[2])
{
  struct diagsrv_srv *srv = (struct diagsrv_srv *) ds->srv;
  struct mrc_fld *f2 = mrc_fld_create(MPI_COMM_SELF);
  mrc_fld_set_param_int_array(f2, "dims", 3, (int[3]) { gdims[0], gdims[1], 1 });
  mrc_fld_set_array(f2, srv->gfld);
  mrc_fld_setup(f2);
  f2->_domain = srv->domain; // FIXME, quite a hack
  return f2;
}

static void
ds_srv_put_gfld_2d(struct diagsrv_one *ds, struct mrc_fld *gfld, char *fld_name,
		   int outtype, float sheet)
{
  mrc_fld_set_comp_name(gfld, 0, fld_name);
  mrc_io_write_field2d(ds->io, 1., gfld, outtype, sheet);
  mrc_fld_destroy(gfld);
}

static void
ds_srv_put_gfld_3d(struct diagsrv_one *ds, struct mrc_fld *gfld)
{
  mrc_fld_write(gfld, ds->io);
  mrc_fld_destroy(gfld);
}

static void
ds_srv_close(struct diagsrv_one *ds)
{
  mrc_io_close(ds->io);
}

static void
ds_srv_write_attr(struct diagsrv_one *ds, const char *path, int type,
		  const char *name, union param_u *pv)
{
  mrc_io_write_attr(ds->io, path, type, name, pv);
}

struct diagsrv_srv_ops ds_srv_ops = {
  .name        = "nocache",
  .create      = ds_srv_create,
  .set_domain  = ds_srv_set_domain,
  .destroy     = ds_srv_destroy,
  .open        = ds_srv_open,
  .get_gfld_2d = ds_srv_get_gfld_2d,
  .get_gfld_3d = ds_srv_get_gfld_3d,
  .put_gfld_2d = ds_srv_put_gfld_2d,
  .put_gfld_3d = ds_srv_put_gfld_3d,
  .write_attr  = ds_srv_write_attr,
  .close       = ds_srv_close,
};

// ----------------------------------------------------------------------
// diagsrv_srv_cache

#define MAX_FIELDS (21)

struct mrc_attrs_entry {
  struct mrc_obj *attrs;
  char *path;
  list_t entry;
};

struct diagsrv_srv_cache_ctx {
  char *obj_names[MAX_FIELDS];
  char *fld_names[MAX_FIELDS];
  int outtypes[MAX_FIELDS];
  float sheets[MAX_FIELDS];
  float *gflds[MAX_FIELDS];
  struct mrc_domain *domain;
  int nr_flds;
  int step;
  float time;
  list_t attrs_list;
};

static void
ds_srv_cache_create(struct diagsrv_one *ds)
{
  struct diagsrv_srv_cache_ctx *srv = calloc(1, sizeof(*srv));
  ds->srv = srv;
}

static void
ds_srv_cache_set_domain(struct diagsrv_one *ds, struct mrc_domain *domain)
{
  struct diagsrv_srv_cache_ctx *srv = ds->srv;
  srv->domain = domain;
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int nglobal = gdims[0] * gdims[1] * gdims[2];
  for (int i = 0; i < MAX_FIELDS; i++) {
    srv->obj_names[i] = NULL;
    srv->fld_names[i] = NULL;
    srv->gflds[i] = malloc(nglobal * sizeof(float));
  }
}

static void
ds_srv_cache_open(struct diagsrv_one *ds, int step, float time)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  srv->nr_flds = 0;
  srv->step = step;
  srv->time = time;
  INIT_LIST_HEAD(&srv->attrs_list);
}

static void
ds_srv_cache_destroy(struct diagsrv_one *ds)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  for (int i = 0; i < MAX_FIELDS; i++) {
    free(srv->obj_names[i]);
    free(srv->fld_names[i]);
    free(srv->gflds[i]);
  }
  free(srv);
}

static struct mrc_fld *
ds_srv_cache_get_gfld_2d(struct diagsrv_one *ds, int gdims[2])
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  free(srv->fld_names[srv->nr_flds]);
  struct mrc_fld *f2 = mrc_fld_create(MPI_COMM_SELF);
  mrc_fld_set_param_int_array(f2, "dims", 3, (int[3]) { gdims[0], gdims[1], 1 });
  mrc_fld_set_array(f2, srv->gflds[srv->nr_flds]);
  mrc_fld_setup(f2);
  f2->_domain = srv->domain; // FIXME, quite a hack
  return f2;
}

static struct mrc_fld *
ds_srv_cache_get_gfld_3d(struct diagsrv_one *ds, int gdims[3])
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  free(srv->obj_names[srv->nr_flds]);
  free(srv->fld_names[srv->nr_flds]);
  struct mrc_fld *fld = mrc_domain_fld_create(srv->domain, SW_0, NULL);
  mrc_fld_set_array(fld, srv->gflds[srv->nr_flds]);
  mrc_fld_setup(fld);
  return fld;
}

static void
ds_srv_cache_put_gfld_2d(struct diagsrv_one *ds, struct mrc_fld *gfld, char *fld_name,
			 int outtype, float sheet)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  srv->fld_names[srv->nr_flds] = strdup(fld_name);
  srv->outtypes[srv->nr_flds] = outtype;
  srv->sheets[srv->nr_flds] = sheet;
  srv->nr_flds++;
  mrc_fld_destroy(gfld);
}

static void
ds_srv_cache_put_gfld_3d(struct diagsrv_one *ds, struct mrc_fld *fld)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  srv->obj_names[srv->nr_flds] = strdup(mrc_fld_name(fld));
  srv->fld_names[srv->nr_flds] = strdup(mrc_fld_comp_name(fld, 0));
  srv->outtypes[srv->nr_flds] = DIAG_TYPE_3D;
  srv->nr_flds++;
  mrc_fld_destroy(fld);
}

static void
ds_srv_cache_write_attr(struct diagsrv_one *ds, const char *path, int type,
			const char *name, union param_u *pv)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;

  struct mrc_attrs_entry *p;
  __list_for_each_entry(p, &srv->attrs_list, entry, struct mrc_attrs_entry) {
    if (strcmp(p->path, path) == 0)
      goto found;
  }

  p = calloc(1, sizeof(*p));
  p->attrs = mrc_obj_create(MPI_COMM_SELF);
  mrc_obj_set_name(p->attrs, path);
  p->path = strdup(path);
  list_add_tail(&p->entry, &srv->attrs_list);

 found:
  mrc_obj_dict_add(p->attrs, type, name, pv);
}

static void
ds_srv_cache_close(struct diagsrv_one *ds)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  int gdims[3];
  mrc_domain_get_global_dims(srv->domain, gdims);

  mrc_io_open(ds->io, "w", srv->step, srv->time);

  for (int i = 0; i < srv->nr_flds; i++) {
    switch (srv->outtypes[i]) {
    case DIAG_TYPE_3D: {
      struct mrc_fld *gfld = mrc_domain_fld_create(srv->domain, SW_0, srv->fld_names[i]);
      mrc_fld_set_array(gfld, srv->gflds[i]);
      mrc_fld_set_name(gfld, srv->obj_names[i]);
      mrc_fld_setup(gfld);
      mrc_fld_write(gfld, ds->io);
      mrc_fld_destroy(gfld);
      break;
    }
    case DIAG_TYPE_2D_X: {
      struct mrc_fld *gfld2 = mrc_fld_create(MPI_COMM_SELF);
      mrc_fld_set_param_int_array(gfld2, "dims", 3, (int[3]) { gdims[1], gdims[2], 1 });
      mrc_fld_set_array(gfld2, srv->gflds[i]);
      mrc_fld_set_comp_name(gfld2, 0, srv->fld_names[i]);
      mrc_fld_setup(gfld2);
      gfld2->_domain = srv->domain; // FIXME, quite a hack
      mrc_io_write_field2d(ds->io, 1., gfld2, DIAG_TYPE_2D_X, srv->sheets[i]);
      mrc_fld_destroy(gfld2);
      break;
    }
    case DIAG_TYPE_2D_Y: {
      struct mrc_fld *gfld2 = mrc_fld_create(MPI_COMM_SELF);
      mrc_fld_set_param_int_array(gfld2, "dims", 3, (int[3]) { gdims[0], gdims[2], 1 });
      mrc_fld_set_array(gfld2, srv->gflds[i]);
      mrc_fld_set_comp_name(gfld2, 0, srv->fld_names[i]);
      mrc_fld_setup(gfld2);
      gfld2->_domain = srv->domain; // FIXME, quite a hack
      mrc_io_write_field2d(ds->io, 1., gfld2, DIAG_TYPE_2D_Y, srv->sheets[i]);
      mrc_fld_destroy(gfld2);
      break;
    }
    case DIAG_TYPE_2D_Z: {
      struct mrc_fld *gfld2 = mrc_fld_create(MPI_COMM_SELF);
      mrc_fld_set_param_int_array(gfld2, "dims", 3, (int[3]) { gdims[0], gdims[1], 1 });
      mrc_fld_set_array(gfld2, srv->gflds[i]);
      mrc_fld_set_comp_name(gfld2, 0, srv->fld_names[i]);
      mrc_fld_setup(gfld2);
      gfld2->_domain = srv->domain; // FIXME, quite a hack
      mrc_io_write_field2d(ds->io, 1., gfld2, DIAG_TYPE_2D_Z, srv->sheets[i]);
      mrc_fld_destroy(gfld2);
      break;
    }
    default:
      assert(0);
    }
  }

  while (!list_empty(&srv->attrs_list)) {
    struct mrc_attrs_entry *p =
      list_entry(srv->attrs_list.next, struct mrc_attrs_entry, entry);
    mrc_obj_write(p->attrs, ds->io);
    mrc_obj_destroy(p->attrs);
    free(p->path);
    list_del(&p->entry);
    free(p);
  }

  mrc_io_close(ds->io);
}

struct diagsrv_srv_ops ds_srv_cache_ops = {
  .name        = "cache",
  .create      = ds_srv_cache_create,
  .destroy     = ds_srv_cache_destroy,
  .set_domain  = ds_srv_cache_set_domain,
  .open        = ds_srv_cache_open,
  .get_gfld_2d = ds_srv_cache_get_gfld_2d,
  .get_gfld_3d = ds_srv_cache_get_gfld_3d,
  .put_gfld_2d = ds_srv_cache_put_gfld_2d,
  .put_gfld_3d = ds_srv_cache_put_gfld_3d,
  .write_attr  = ds_srv_cache_write_attr,
  .close       = ds_srv_cache_close,
};

// ----------------------------------------------------------------------
// diagsrv_one helpers

static void
add_to_field_2d(struct mrc_fld *g, struct mrc_fld *l, int ib[2])
{
  struct mrc_fld *_g = mrc_fld_get_as(g, "float");
  struct mrc_fld *_l = mrc_fld_get_as(l, "float");
  for (int iy = 0; iy < l->_dims.vals[1]; iy++) {
    for (int ix = 0; ix < l->_dims.vals[0]; ix++) {
      MRC_F2(_g,0, ix+ib[0],iy+ib[1]) = MRC_F2(_l,0, ix,iy);
    }
  }
  mrc_fld_put_as(_g, g);
  mrc_fld_put_as(_l, l);
}

static void
add_to_field_3d(struct mrc_fld *g_, struct mrc_fld *l_, int ib[3])
{
  struct mrc_fld *g = mrc_fld_get_as(g_, "float");
  struct mrc_fld *l = mrc_fld_get_as(l_, "float");
  mrc_fld_foreach(l, ix,iy,iz, 0, 0) {
    MRC_F3(g,0, ix+ib[0],iy+ib[1],iz+ib[2]) = MRC_F3(l,0, ix,iy,iz);
  } mrc_fld_foreach_end;
  mrc_fld_put_as(g, g_);
  mrc_fld_put_as(l, l_);
}

// ----------------------------------------------------------------------
// diagsrv_one

static struct diagsrv_srv_ops *ds_srvs[] = {
  &ds_srv_ops,
  &ds_srv_cache_ops,
  NULL,
};

static struct diagsrv_srv_ops *
find_ds_srv(const char *ds_srv)
{
#if 0
  fprintf(stderr,"Available ds_srvs:\n");
  for (int i = 0; ds_srvs[i]; i++){
    fprintf(stderr,"ds_srvs[i]:%s\n",ds_srvs[i]->name);
  }
#endif

  for (int i = 0; ds_srvs[i]; i++) {
    if (strcmp(ds_srv, ds_srvs[i]->name) == 0) {
      return ds_srvs[i];
    }
  }
  fprintf(stderr, "ERROR: unknown diagsrv_srv (ds_srv) '%s'\n", ds_srv);
  abort();
}


static struct mrc_domain *
diagsrv_recv_domain_info(int nr_procs, int ldims[3])
{
  for (int d = 0; d < 3; d++) {
    ldims[d] = 0;
  }

  struct mrc_domain *domain = NULL;
  struct mrc_crds *crds = NULL;
  int iw[9], *off = iw, *_ldims = iw + 3, *gdims = iw + 6;
  for (int rank = 0; rank < nr_procs; rank++) {
    MPI_Recv(iw, 9, MPI_INT, rank, ID_DIAGS_CMD_DOMAIN_INFO, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (rank == 0) {
      domain = mrc_domain_create(MPI_COMM_SELF);
      mrc_domain_set_type(domain, "simple");
      mrc_domain_set_param_int3(domain, "m", gdims);
      crds = mrc_domain_get_crds(domain);
      mrc_crds_set_type(crds, "rectilinear");
      mrc_domain_setup(domain);
    }
    for (int d = 0; d < 3; d++) {
      // find max local domain
      if (ldims[d] < _ldims[d]) {
        ldims[d] = _ldims[d];
      }
      // FIXME, this is a bit stupid, but by far the easiest way without assuming
      // internal mpi_domain knowledge -- we repeatedly receive the same pieces
      // of the global coord array, but it's only done once, anyway
      float *buf = calloc(_ldims[d], sizeof(*buf));
      MPI_Recv(buf, _ldims[d], MPI_FLOAT, rank, ID_DIAGS_CMD_CRDX + d, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
      for (int i = 0; i < _ldims[d]; i++) {
        MRC_CRD(crds, d, i + off[d]) = buf[i];
      }
      free(buf);
    }
  }
  return domain;
}

void
mrc_io_server(const char *ds_format, const char *ds_srv, int nr_procs)
{
  struct diagsrv_params {
    char *format;
    char *server;
  };

#define VAR(x) (void *)offsetof(struct diagsrv_params, x)
static struct param diagsrv_params_descr[] = {
  { "diagsrv_format"     , VAR(format)          , PARAM_STRING(NULL)      },
  { "diagsrv_server"     , VAR(server)          , PARAM_STRING(NULL)      },
  {},
};
#undef VAR

  struct diagsrv_params par = {
    .format = (char *) ds_format,
    .server = (char *) ds_srv,
  };
  mrc_params_parse_nodefault(&par, diagsrv_params_descr, "diagsrv",
			     MPI_COMM_SELF);
  mrc_params_print(&par, diagsrv_params_descr, "diagsrv", MPI_COMM_SELF);

  struct diagsrv_srv_ops *srv_ops = find_ds_srv(par.server);
  struct diagsrv_one ds = {
    .srv_ops    = srv_ops,
  };
  INIT_LIST_HEAD(&ds.mrc_io_list);

  int respond_to_rank = -1;
  int gdims[3];
  float *w2 = NULL;

  srv_ops->create(&ds);

  // loop waiting for data to write
  for (;;) {
    int icmd[2];
    MPI_Status status;
    MPI_Status stat;
    MPI_Recv(icmd, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    if (icmd[0] == DIAG_CMD_SHUTDOWN) {
      if (icmd[1]) {
        respond_to_rank = stat.MPI_SOURCE;
      }
      break;
    }

    if (icmd[0] == DIAG_CMD_CREATE) {
      char s[256] = {};
      MPI_Recv(s, 255, MPI_CHAR, 0, ID_DIAGS_CREATE_OUTDIR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      char *outdir = strdup(s);
      MPI_Recv(s, 255, MPI_CHAR, 0, ID_DIAGS_CREATE_BASENAME, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      char *basename = strdup(s);

      assert(!find_diagsrv_io(&ds, basename));
      create_diagsrv_io(&ds, par.format, outdir, basename);
      continue;
    }

    assert(icmd[0] == DIAG_CMD_OPENFILE);
    int step = icmd[1];

    int outtag = stat.MPI_TAG;

    char s[256] = {};
    MPI_Recv(s, 255, MPI_CHAR, 0, ID_DIAGS_BASENAME, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    struct mrc_io_entry *io_entry = find_diagsrv_io(&ds, s);
    assert(io_entry);
    ds.io = io_entry->io;

    float time;
    MPI_Recv(&time, 1, MPI_FLOAT, 0, ID_DIAGS_TIME, MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);

    switch (outtag) {
    case ID_DIAGS_CMD_OPEN:
      srv_ops->open(&ds, step, time);
      break;
    default:
      assert(0);
    }

    for (;;) { //waiting for field to write.

      char fld_name80[80];
      MPI_Recv(fld_name80, 80, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
	       &status);
      if (status.MPI_TAG == ID_DIAGS_CMD_WRITE_ATTR) {
	char *path = fld_name80;
	int type;
	MPI_Recv(&type, 1, MPI_INT, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	char name[80];
	MPI_Recv(name, 80, MPI_CHAR, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	union param_u val;
	switch (type) {
	case PT_BOOL:
	case PT_INT:
	case PT_SELECT:
	case MRC_VAR_INT:
	case MRC_VAR_BOOL:
	  MPI_Recv(&val.u_int, 1, MPI_INT, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  break;
	case PT_FLOAT:
	  MPI_Recv(&val.u_float, 1, MPI_FLOAT, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  break;
	case PT_DOUBLE:
	case MRC_VAR_DOUBLE:
	  MPI_Recv(&val.u_double, 1, MPI_DOUBLE, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  break;
	case PT_STRING: ;
	  char str[100];
	  MPI_Recv(str, 100, MPI_CHAR, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  val.u_string = strdup(str);
	  break;
	case PT_INT3:
	  MPI_Recv(val.u_int3, 3, MPI_INT, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  break;
	case PT_FLOAT3:
	  MPI_Recv(val.u_float3, 3, MPI_FLOAT, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  break;
	case PT_DOUBLE3:
	case MRC_VAR_DOUBLE3:
	  MPI_Recv(val.u_double3, 3, MPI_DOUBLE, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  break;
	case PT_INT_ARRAY:
	  MPI_Recv(&val.u_int_array.nr_vals, 1, MPI_INT, 0, ID_DIAGS_CMD_WRITE_ATTR, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  val.u_int_array.vals = calloc(val.u_int_array.nr_vals, sizeof(int));
	  MPI_Recv(val.u_int_array.vals, val.u_int_array.nr_vals, MPI_INT, 0, ID_DIAGS_CMD_WRITE_ATTR,
		   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  break;
	default:
	  assert(0);
	}
	srv_ops->write_attr(&ds, path, type, name, &val);
	continue;
      } else if (status.MPI_TAG == ID_DIAGS_CMD_CREATE) {
	assert (!io_entry->domain);
	io_entry->domain = diagsrv_recv_domain_info(nr_procs, io_entry->ldims);
	srv_ops->set_domain(&ds, io_entry->domain);
	mrc_domain_get_global_dims(io_entry->domain, gdims);
	int *ldims = io_entry->ldims;
	w2 = malloc(ldims[0] * ldims[1] * ldims[2] * sizeof(float));
	continue;
      };
      assert(status.MPI_TAG == ID_DIAGS_FLDNAME);
      assert(io_entry->domain);

      if (!fld_name80[0])
	break;

      char obj_name80[80];
      strcpy(obj_name80, fld_name80);

      MPI_Recv(fld_name80, 80, MPI_CHAR, 0, ID_DIAGS_FLDNAME, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      int outtype;
      MPI_Recv(&outtype, 1, MPI_INT, 0, ID_DIAGS_CMD_WRITE, MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);

      if (outtype != DIAG_TYPE_3D) {
	float sheet;
	MPI_Recv(&sheet, 1, MPI_FLOAT, 0, ID_DIAGS_CMD_WRITE, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);

	int i0 = -1, i1 = -1;
	switch (outtype) {
	case DIAG_TYPE_2D_X: i0 = 1; i1 = 2; break;
	case DIAG_TYPE_2D_Y: i0 = 0; i1 = 2; break;
	case DIAG_TYPE_2D_Z: i0 = 0; i1 = 1; break;
	default:
	  assert(0);
	}

	struct mrc_fld *gfld2 = srv_ops->get_gfld_2d(&ds, (int [2]) { gdims[i0], gdims[i1] });

	for (int k = 0; k < nr_procs; k++) {
	  int iw[6], *off = iw, *dims = iw + 3; // off, then dims
	  MPI_Recv(iw, 6, MPI_INT, k, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  // receive data and add to field
	  if (iw[0] > -1) {
	    struct mrc_fld *lfld2 = mrc_fld_create(MPI_COMM_SELF);
	    mrc_fld_set_param_int_array(lfld2, "dims", 3, (int[3]) { dims[i0], dims[i1], 1 });
	    mrc_fld_set_array(lfld2, w2);
	    mrc_fld_setup(lfld2);
	    struct mrc_fld *_lfld2 = mrc_fld_get_as(lfld2, "float");
	    MPI_Recv(_lfld2->_nd->arr, mrc_fld_len(lfld2), MPI_FLOAT, k, ID_DIAGS_2DDATA, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    mrc_fld_put_as(_lfld2, lfld2);
	    add_to_field_2d(gfld2, lfld2, (int [2]) { off[i0], off[i1] });
	    mrc_fld_destroy(lfld2);
	  }
	}

	srv_ops->put_gfld_2d(&ds, gfld2, fld_name80, outtype, sheet);
      } else {
	struct mrc_fld *gfld3 = srv_ops->get_gfld_3d(&ds, gdims);

	for (int k = 0; k < nr_procs; k++) {
	  int iw[6], *off = iw, *dims = iw + 3; // off, then dims
	  MPI_Recv(iw, 6, MPI_INT, k, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  // receive data and add to field
	  if (iw[0] > -1) {
	    struct mrc_fld *lfld3 = mrc_fld_create(MPI_COMM_SELF);
	    mrc_fld_set_param_int_array(lfld3, "dims", 4,
				       (int[4]) { dims[0], dims[1], dims[2], 1 });
	    mrc_fld_set_array(lfld3, w2);
	    mrc_fld_setup(lfld3);
	    struct mrc_fld *_lfld3 = mrc_fld_get_as(lfld3, "float");
	    MPI_Recv(lfld3->_nd->arr, mrc_fld_len(lfld3), MPI_FLOAT, k, ID_DIAGS_DATA, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    mrc_fld_put_as(_lfld3, lfld3);
	    add_to_field_3d(gfld3, lfld3, off);
	    mrc_fld_destroy(lfld3);
	  }
	}

	mrc_fld_set_name(gfld3, obj_name80);
	mrc_fld_set_comp_name(gfld3, 0, fld_name80);
	srv_ops->put_gfld_3d(&ds, gfld3);
      }
    }
    srv_ops->close(&ds);
    ds.io = NULL;
  }  //for (;;) //loop waiting for data to write...
  free(w2);

  while (!list_empty(&ds.mrc_io_list)) {
    struct mrc_io_entry *p = list_entry(ds.mrc_io_list.next, struct mrc_io_entry, entry);
    mrc_io_destroy(p->io);
    mrc_domain_destroy(p->domain);
    list_del(&p->entry);
    free(p);
  }

  srv_ops->destroy(&ds);

  if (respond_to_rank > -1) {
    int response = DIAG_RESPONSE_SHUTDOWN_COMPLETE;
    assert(respond_to_rank == 0);
    MPI_Send(&response, 1, MPI_INT, respond_to_rank, ID_DIAGS_CMD_RESPONSE,
	     MPI_COMM_WORLD);
  }
}
