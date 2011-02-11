
#include "mrc_diag_private.h"
#include <mrc_params.h>

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
  ID_DIAGS_CMD_OUTDIR,
  ID_DIAGS_CMD_RUN,
  ID_DIAGS_CMD_CREATE,
  ID_DIAGS_CMD_CRDX,
  ID_DIAGS_CMD_CRDY,
  ID_DIAGS_CMD_CRDZ,
  ID_DIAGS_CMD_CREATE_3D,
  ID_DIAGS_CMD_CREATE_X,
  ID_DIAGS_CMD_CREATE_Y,
  ID_DIAGS_CMD_CREATE_Z,
  ID_DIAGS_TIME,
  ID_DIAGS_SHEET,
  ID_DIAGS_TIMESTR,
  ID_DIAGS_FLDNAME,
  ID_DIAGS_SUBDOMAIN,
  ID_DIAGS_DATA,
  ID_DIAGS_2DDATA,
};

enum {
  DIAG_CMD_OPENFILE,
  DIAG_CMD_SHUTDOWN,
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
diagc_combined_send_domain_info(struct diag_format *format, struct mrc_domain *domain)
{
  if (format->diagc_domain_info_sent)
    return;

  struct diagc_combined_params *par = format->obj.subctx;
    
  int iw[9], *off = iw, *ldims = iw + 3, *gdims = iw + 6;
  mrc_domain_get_local_offset_dims(domain, off, ldims);
  mrc_domain_get_global_dims(domain, gdims);
  MPI_Send(iw, 9, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD_CREATE, MPI_COMM_WORLD);

  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  for (int d = 0; d < 3; d++) {
    MPI_Send(crds->crd[d] + 2, ldims[d], MPI_FLOAT, par->rank_diagsrv,
	     ID_DIAGS_CMD_CRDX + d, MPI_COMM_WORLD); // FIXME, hardcoded SW
  }

  format->diagc_domain_info_sent = true;
}

static void
diagc_combined_setup(struct mrc_obj *obj)
{
  struct diag_format *format = to_diag_format(obj);
  struct diag_info *info = &format->diag_info;
  struct diagc_combined_params *par = format->obj.subctx;

  if (format->rank == 0) {
    MPI_Send((char *)info->outdir, strlen(info->outdir), MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_CMD_OUTDIR, MPI_COMM_WORLD);
    MPI_Send((char *)info->run, strlen(info->run), MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_CMD_RUN, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_open
//
// opens a new output file (index itime)
// afterwards, diag_write_field() may be called repeatedly and
// must be concluded with a call to diag_close()

static void
diagc_combined_open(struct diag_format *format, float sheet, int outtype, int step)
{
  struct diagc_combined_params *par = format->obj.subctx;

  if (format->rank == 0) {
    int icmd[2] = { DIAG_CMD_OPENFILE, step };
    int tag;

    switch (outtype) {
    case DIAG_TYPE_3D  : tag = ID_DIAGS_CMD_CREATE_3D; break;
    case DIAG_TYPE_2D_X: tag = ID_DIAGS_CMD_CREATE_X; break;
    case DIAG_TYPE_2D_Y: tag = ID_DIAGS_CMD_CREATE_Y; break;
    case DIAG_TYPE_2D_Z: tag = ID_DIAGS_CMD_CREATE_Z; break;
    default: assert(0);
    }

    MPI_Send(icmd, 2, MPI_INT, par->rank_diagsrv, tag, MPI_COMM_WORLD);
    MPI_Send(&format->time, 1, MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_TIME, MPI_COMM_WORLD);
    MPI_Send(format->time_str, strlen(format->time_str), MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_TIMESTR, MPI_COMM_WORLD);
    MPI_Send(&sheet, 1, MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_SHEET, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_close

static void
diagc_combined_close(struct diag_format *format)
{
  struct diagc_combined_params *par = format->obj.subctx;

  if (format->rank == 0) {
    char str[] = "";
    MPI_Send(str, 1, MPI_CHAR, par->rank_diagsrv, ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_destroy
//
// shuts down the diag server process

static void
diagc_combined_destroy(struct mrc_obj *obj)
{
  struct diag_format *format = (struct diag_format *) obj;
  struct diagc_combined_params *par = format->obj.subctx;

  if (format->rank == 0) {
    int icmd[1] = { DIAG_CMD_SHUTDOWN };
    MPI_Send(icmd, 1, MPI_INT, par->rank_diagsrv, ID_DIAGS_CMD_CREATE_3D, MPI_COMM_WORLD);
  }
}

// ----------------------------------------------------------------------
// diagc_combined_write_field

static void
copy_and_scale(float *buf, struct mrc_f3 *fld, int m, float scale)
{
  int i = 0;
  for (int iz = 2; iz < fld->im[2] - 2; iz++) {
    for (int iy = 2; iy < fld->im[1] - 2; iy++) {
      for (int ix = 2; ix < fld->im[0] - 2; ix++) {
	buf[i++] = scale * MRC_F3(fld, m, ix,iy,iz);
      }
    }
  }
}

static void
diagc_combined_write_field(struct diag_format *format, float scale,
			   struct mrc_f3 *fld, int m)
{
  struct diagc_combined_params *par = format->obj.subctx;

  diagc_combined_send_domain_info(format, fld->domain);

  int ldims[3];
  mrc_domain_get_local_offset_dims(fld->domain, NULL, ldims);
  assert(ldims[0] == fld->im[0]-4 && ldims[1] == fld->im[1]-4 && ldims[2] == fld->im[2]-4);
  int nout = ldims[0] * ldims[1] * ldims[2];
  float *buf = calloc(sizeof(float), nout);

  if (format->rank == 0) {
    MPI_Send(fld->name[m], strlen(fld->name[m]), MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
  }

  copy_and_scale(buf, fld, m, scale);

  int iw[6], *off = iw, *dims = iw + 3; // off, then dims
  mrc_domain_get_local_offset_dims(fld->domain, off, dims);

  MPI_Send(iw, 6, MPI_INT, par->rank_diagsrv, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD);
  MPI_Send(buf, nout, MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_DATA, MPI_COMM_WORLD);

  free(buf);
}

static void
diagc_combined_write_field2d(struct diag_format *format, float scale, struct mrc_f2 *fld,
			     const char *fld_name, int outtype, float sheet)
{
  struct diagc_combined_params *par = format->obj.subctx;

  diagc_combined_send_domain_info(format, fld->domain);

  assert(outtype >= DIAG_TYPE_2D_X && outtype <= DIAG_TYPE_2D_Z);
  int dim = outtype - DIAG_TYPE_2D_X;

  if (format->rank == 0) {
    MPI_Send((char *)fld_name, strlen(fld_name), MPI_CHAR, par->rank_diagsrv,
	     ID_DIAGS_FLDNAME, MPI_COMM_WORLD);
  }

  int iw[6] = { -1, };
  if (fld->arr) { // part of the slice?
    int *off = iw, *dims = iw + 3; // off, then dims
    mrc_domain_get_local_offset_dims(fld->domain, off, dims);
    off[dim] = 0;
    dims[dim] = 1;

    MPI_Send(iw, 6, MPI_INT, par->rank_diagsrv, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD);
    MPI_Send(fld->arr, fld->len, MPI_FLOAT, par->rank_diagsrv, ID_DIAGS_2DDATA, MPI_COMM_WORLD);
  } else {
    MPI_Send(iw, 6, MPI_INT, par->rank_diagsrv, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD);
  }
}

static struct diag_format_ops ds_combined_ops = {
  .name          = "combined",
  .size          = sizeof(struct diagc_combined_params),
  .param_descr   = diagc_combined_params_descr,
  .setup         = diagc_combined_setup,
  .destroy       = diagc_combined_destroy,
  .open          = diagc_combined_open,
  .write_field   = diagc_combined_write_field,
  .write_field2d = diagc_combined_write_field2d,
  .close         = diagc_combined_close,
};

// ======================================================================
// diag server one

// generic diagnostics server running on one node.
// a particular implementation needs to take care of opening files/
// writing/closing them,
// the generic part here (in conjunction with the client side above)
// takes care of communicating / assembling the fields

struct diagsrv_one {
  struct diag_format *format;
  struct diagsrv_srv_ops *srv_ops;
  void *srv;
};

struct diagsrv_srv_ops {
  const char *name;
  void  (*create)(struct diagsrv_one *ds, struct mrc_domain *domain);
  void  (*destroy)(struct diagsrv_one *ds);
  void  (*open)(struct diagsrv_one *ds, float sheet, int outtype, int step,
		float time, const char *time_str);
  void  (*get_gfld_2d)(struct diagsrv_one *ds, struct mrc_f2 *f2, int dims[2]);
  void  (*get_gfld_3d)(struct diagsrv_one *ds, struct mrc_f3 *f3, int dims[3]);
  void  (*put_gfld_2d)(struct diagsrv_one *ds, struct mrc_f2 *f2, char *fld_name, 
		       int outtype, float sheet);
  void  (*put_gfld_3d)(struct diagsrv_one *ds, struct mrc_f3 *f3);
  void  (*close)(struct diagsrv_one *ds);
};

// ----------------------------------------------------------------------
// diagsrv_one

struct diagsrv_srv {
  float *gfld;
  struct mrc_domain *domain;
};

static void
ds_srv_create(struct diagsrv_one *ds, struct mrc_domain *domain)
{
  struct diagsrv_srv *srv = malloc(sizeof(*srv));
  srv->domain = domain;
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int nglobal = gdims[0] * gdims[1] * gdims[2];
  srv->gfld = malloc(nglobal * sizeof(float));
  ds->srv = srv;
}

static void
ds_srv_destroy(struct diagsrv_one *ds)
{
  struct diagsrv_srv *srv = (struct diagsrv_srv *) ds->srv;
  free(srv->gfld);
  free(srv);
}

static void
ds_srv_open(struct diagsrv_one *ds, float sheet, int outtype, int step,
	    float time, const char *time_str)
{
  diag_format_open(ds->format, sheet, outtype, step, time, time_str);
}

static void
ds_srv_get_gfld_3d(struct diagsrv_one *ds, struct mrc_f3 *f3, int gdims[3])
{
  struct diagsrv_srv *srv = (struct diagsrv_srv *) ds->srv;
  mrc_domain_f3_alloc_with_array(srv->domain, f3, 1, SW_0, srv->gfld);
}

static void
ds_srv_get_gfld_2d(struct diagsrv_one *ds, struct mrc_f2 *f2, int gdims[2])
{
  struct diagsrv_srv *srv = (struct diagsrv_srv *) ds->srv;
  mrc_f2_alloc_with_array(f2, NULL, gdims, 1, srv->gfld);
  f2->domain = srv->domain; // FIXME, quite a hack
}

static void
ds_srv_put_gfld_2d(struct diagsrv_one *ds, struct mrc_f2 *gfld, char *fld_name,
		   int outtype, float sheet)
{
  diag_format_write_field2d(ds->format, 1., gfld, fld_name, outtype, sheet);
}

static void
ds_srv_put_gfld_3d(struct diagsrv_one *ds, struct mrc_f3 *gfld)
{
  diag_format_write_field(ds->format, 1., gfld, 0);
  mrc_f3_free(gfld);
}

static void
ds_srv_close(struct diagsrv_one *ds)
{
  diag_format_close(ds->format);
}

struct diagsrv_srv_ops ds_srv_ops = {
  .name        = "nocache",
  .create      = ds_srv_create,
  .destroy     = ds_srv_destroy,
  .open        = ds_srv_open,
  .get_gfld_2d = ds_srv_get_gfld_2d,
  .get_gfld_3d = ds_srv_get_gfld_3d,
  .put_gfld_2d = ds_srv_put_gfld_2d,
  .put_gfld_3d = ds_srv_put_gfld_3d,
  .close       = ds_srv_close,
};

// ----------------------------------------------------------------------
// diagsrv_srv_cache

#define MAX_FIELDS (21)

struct diagsrv_srv_cache_ctx {
  char *fld_names[MAX_FIELDS];
  float sheets[MAX_FIELDS];
  float *gflds[MAX_FIELDS];
  struct mrc_domain *domain;
  int nr_flds;
  float sheet;
  int outtype;
  int step;
  float time;
  char *time_str;
};

static void
ds_srv_cache_create(struct diagsrv_one *ds, struct mrc_domain *domain)
{
  struct diagsrv_srv_cache_ctx *srv = malloc(sizeof(*srv));
  srv->domain = domain;
  int gdims[3];
  mrc_domain_get_global_dims(domain, gdims);
  int nglobal = gdims[0] * gdims[1] * gdims[2];
  for (int i = 0; i < MAX_FIELDS; i++) {
    srv->fld_names[i] = NULL;
    srv->gflds[i] = malloc(nglobal * sizeof(float));
  }
  ds->srv = srv;
}

static void
ds_srv_cache_open(struct diagsrv_one *ds, float sheet, int outtype, int step,
		  float time, const char *time_str)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  srv->nr_flds = 0;
  srv->sheet = sheet;
  srv->outtype = outtype;
  srv->step = step;
  srv->time = time;
  srv->time_str = strdup(time_str);
}

static void
ds_srv_cache_destroy(struct diagsrv_one *ds)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  for (int i = 0; i < MAX_FIELDS; i++) {
    free(srv->fld_names[i]);
    free(srv->gflds[i]);
  }
  free(srv);
}

static void
ds_srv_cache_get_gfld_2d(struct diagsrv_one *ds, struct mrc_f2 *f2, int gdims[2])
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  free(srv->fld_names[srv->nr_flds]);
  mrc_f2_alloc_with_array(f2, NULL, gdims, 1, srv->gflds[srv->nr_flds]);
  f2->domain = srv->domain;
}

static void
ds_srv_cache_get_gfld_3d(struct diagsrv_one *ds, struct mrc_f3 *f3, int gdims[3])
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  free(srv->fld_names[srv->nr_flds]);
  mrc_domain_f3_alloc_with_array(srv->domain, f3, 1, SW_0, srv->gflds[srv->nr_flds]);
}

static void
ds_srv_cache_put_gfld_2d(struct diagsrv_one *ds, struct mrc_f2 *gfld, char *fld_name,
			 int outtype, float sheet)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  srv->fld_names[srv->nr_flds] = strdup(fld_name);
  srv->sheets[srv->nr_flds] = sheet;
  srv->nr_flds++;
}

static void
ds_srv_cache_put_gfld_3d(struct diagsrv_one *ds, struct mrc_f3 *f3)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  assert(srv->nr_flds < MAX_FIELDS);
  srv->fld_names[srv->nr_flds] = strdup(f3->name[0]);
  srv->nr_flds++;
  mrc_f3_free(f3);
}

static void
ds_srv_cache_close(struct diagsrv_one *ds)
{
  struct diagsrv_srv_cache_ctx *srv = (struct diagsrv_srv_cache_ctx *) ds->srv;
  int gdims[3];
  mrc_domain_get_global_dims(srv->domain, gdims);

  diag_format_open(ds->format, srv->sheet, srv->outtype, srv->step, srv->time, srv->time_str);

  for (int i = 0; i < srv->nr_flds; i++) {
    switch (srv->outtype) {
    case DIAG_TYPE_3D: {
      struct mrc_f3 gfld;
      mrc_domain_f3_alloc_with_array(srv->domain, &gfld, 1, SW_0, srv->gflds[i]);
      gfld.name[0] = strdup(srv->fld_names[i]);
      diag_format_write_field(ds->format, 1., &gfld, 0);
      mrc_f3_free(&gfld);
      break;
    }
    case DIAG_TYPE_2D_X: {
      struct mrc_f2 gfld2;
      mrc_f2_alloc_with_array(&gfld2, NULL, (int [2]) { gdims[1], gdims[2] }, 1, srv->gflds[i]);
      gfld2.domain = srv->domain;
      diag_format_write_field2d(ds->format, 1., &gfld2, srv->fld_names[i], DIAG_TYPE_2D_X, srv->sheets[i]);
      mrc_f2_free(&gfld2);
      break;
    }
    case DIAG_TYPE_2D_Y: {
      struct mrc_f2 gfld2;
      mrc_f2_alloc_with_array(&gfld2, NULL, (int [2]) { gdims[0], gdims[2] }, 1, srv->gflds[i]);
      gfld2.domain = srv->domain;
      diag_format_write_field2d(ds->format, 1., &gfld2, srv->fld_names[i], DIAG_TYPE_2D_Y, srv->sheets[i]);
      mrc_f2_free(&gfld2);
      break;
    }
    case DIAG_TYPE_2D_Z: {
      struct mrc_f2 gfld2;
      mrc_f2_alloc_with_array(&gfld2, NULL, (int [2]) { gdims[0], gdims[1] }, 1, srv->gflds[i]);
      gfld2.domain = srv->domain;
      diag_format_write_field2d(ds->format, 1., &gfld2, srv->fld_names[i], DIAG_TYPE_2D_Z, srv->sheets[i]);
      mrc_f2_free(&gfld2);
      break;
    }
    default:
      assert(0);
    }
  }

  diag_format_close(ds->format);
  free(srv->time_str);
}

struct diagsrv_srv_ops ds_srv_cache_ops = {
  .name        = "cache",
  .create      = ds_srv_cache_create,
  .destroy     = ds_srv_cache_destroy,
  .open        = ds_srv_cache_open,
  .get_gfld_2d = ds_srv_cache_get_gfld_2d,
  .get_gfld_3d = ds_srv_cache_get_gfld_3d,
  .put_gfld_2d = ds_srv_cache_put_gfld_2d,
  .put_gfld_3d = ds_srv_cache_put_gfld_3d,
  .close       = ds_srv_cache_close,
};

// ----------------------------------------------------------------------
// diagsrv_one helpers

static void
add_to_field_2d(struct mrc_f2 *g, struct mrc_f2 *l, int ib[2])
{
  for (int iy = 0; iy < l->im[1]; iy++) {
    for (int ix = 0; ix < l->im[0]; ix++) {
      MRC_F2(g,0, ix+ib[0],iy+ib[1]) = MRC_F2(l,0, ix,iy);
    }
  }
}

static void
add_to_field_3d(struct mrc_f3 *g, struct mrc_f3 *l, int ib[2])
{
  for (int iz = 0; iz < l->im[2]; iz++) {
    for (int iy = 0; iy < l->im[1]; iy++) {
      for (int ix = 0; ix < l->im[0]; ix++) {
	MRC_F3(g,0, ix+ib[0],iy+ib[1],iz+ib[2]) = MRC_F3(l,0, ix,iy,iz);
      }
    }
  }
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
    MPI_Recv(iw, 9, MPI_INT, rank, ID_DIAGS_CMD_CREATE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (rank == 0) {
      struct mrc_domain_simple_params par = {
	.ldims    = { gdims[0], gdims[1], gdims[2] },
	.nr_procs = { 1, 1, 1 },
      };
      domain = mrc_domain_create(MPI_COMM_SELF, "simple");
      mrc_domain_simple_set_params(domain, &par);
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
	crds->crd[d][i + off[d]] = buf[i];
      }
      free(buf);
    }
  }
  return domain;
}

void
diagsrv_one(const char *ds_format, const char *ds_srv, int nr_procs)
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
  char *outdir, *run;

  // initial handshake

  char s[256] = {};
  MPI_Recv(s, 255, MPI_CHAR, 0, ID_DIAGS_CMD_OUTDIR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  outdir = strdup(s);
  MPI_Recv(s, 255, MPI_CHAR, 0, ID_DIAGS_CMD_RUN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  run = strdup(s);

  struct diagsrv_one ds = {
    .srv_ops    = srv_ops,
  };
  ds.format = diag_format_create(MPI_COMM_SELF, par.format);
  diag_format_set_param_string(ds.format, "outdir", outdir);
  diag_format_set_param_string(ds.format, "run", run);
  diag_format_view(ds.format);
  diag_format_setup(ds.format);

  float *w2 = NULL;
  struct mrc_domain *domain = NULL;
  int gdims[3];

  // loop waiting for data to write
  for (;;) {
    int icmd[2];
    MPI_Status status;
    MPI_Status stat;
    MPI_Recv(icmd, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    if (icmd[0] == DIAG_CMD_SHUTDOWN) {
      break;
    }
    int outtag = stat.MPI_TAG;

    assert(icmd[0] == DIAG_CMD_OPENFILE);
    int step = icmd[1];

    float time;
    MPI_Recv(&time, 1, MPI_FLOAT, 0, ID_DIAGS_TIME, MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);

    int len;
    char time_str80[81];
    MPI_Recv(time_str80, 80, MPI_CHAR, 0, ID_DIAGS_TIMESTR, MPI_COMM_WORLD,
	     &status);
    MPI_Get_count(&status, MPI_CHAR, &len);
    time_str80[len] = 0;
    
    float sheet;
    MPI_Recv( &sheet, 1, MPI_FLOAT, 0, ID_DIAGS_SHEET, MPI_COMM_WORLD,
	      MPI_STATUS_IGNORE);

    int i0 = -1, i1 = -1;

    if (!w2) {
      int ldims[3];
      domain = diagsrv_recv_domain_info(nr_procs, ldims);
      mrc_domain_get_global_dims(domain, gdims);
      srv_ops->create(&ds, domain);
      w2 = malloc(ldims[0] * ldims[1] * ldims[2] * sizeof(float));
    }

    switch (outtag) {
    case ID_DIAGS_CMD_CREATE_3D:
      srv_ops->open(&ds, 0., DIAG_TYPE_3D, step, time, time_str80);
      break;
    case ID_DIAGS_CMD_CREATE_X: // Assemble a constant X sheet
      srv_ops->open(&ds, sheet, DIAG_TYPE_2D_X, step, time, time_str80);
      i0 = 1; i1 = 2;
      break;
    case ID_DIAGS_CMD_CREATE_Y: // Assemble a constant Y sheet
      srv_ops->open(&ds, sheet, DIAG_TYPE_2D_Y, step, time, time_str80);	
      i0 = 0; i1 = 2;
      break;
    case ID_DIAGS_CMD_CREATE_Z: // Assemble a constant Z sheet
      srv_ops->open(&ds, sheet, DIAG_TYPE_2D_Z, step, time, time_str80);
      i0 = 0; i1 = 1;
      break;
    default:
      assert(0);
    }

    for (;;) { //waiting for field to write.
  
      char fld_name80[81];
      MPI_Recv(fld_name80, 80, MPI_CHAR, 0, ID_DIAGS_FLDNAME, MPI_COMM_WORLD,
	       &status);
      MPI_Get_count(&status, MPI_CHAR, &len);
      fld_name80[len] = 0;
      
      if (!fld_name80[0])
	break;
      
      if (outtag != ID_DIAGS_CMD_CREATE_3D) {
	struct mrc_f2 gfld2, lfld2;
	srv_ops->get_gfld_2d(&ds, &gfld2, (int [2]) { gdims[i0], gdims[i1] });

	for (int k = 0; k < nr_procs; k++) { 
	  int iw[6], *off = iw, *dims = iw + 3; // off, then dims
	  MPI_Recv(iw, 6, MPI_INT, k, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  
	  // receive data and add to field
	  if (iw[0] > -1) {
	    mrc_f2_alloc_with_array(&lfld2, NULL, (int [2]) { dims[i0], dims[i1] }, 1, w2);
	    MPI_Recv(lfld2.arr, lfld2.len, MPI_FLOAT, k, ID_DIAGS_2DDATA, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    add_to_field_2d(&gfld2, &lfld2, (int [2]) { off[i0], off[i1] });
	    mrc_f2_free(&lfld2);
	  }
	}
	
	srv_ops->put_gfld_2d(&ds, &gfld2, fld_name80,
			     outtag - ID_DIAGS_CMD_CREATE_X + DIAG_TYPE_2D_X, sheet);
      } else {
	struct mrc_f3 gfld3, lfld3;
	srv_ops->get_gfld_3d(&ds, &gfld3, gdims);

	for (int k = 0; k < nr_procs; k++) { 
	  int iw[6], *off = iw, *dims = iw + 3; // off, then dims
	  MPI_Recv(iw, 6, MPI_INT, k, ID_DIAGS_SUBDOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  
	  // receive data and add to field
	  if (iw[0] > -1) {
	    mrc_f3_alloc_with_array(&lfld3, NULL, dims, 1, w2);
	    MPI_Recv(lfld3.arr, lfld3.len, MPI_FLOAT, k, ID_DIAGS_DATA, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    add_to_field_3d(&gfld3, &lfld3, off);
	    mrc_f3_free(&lfld3);
	  }
	}

	gfld3.name[0] = strdup(fld_name80);
	srv_ops->put_gfld_3d(&ds, &gfld3);
      }
    }  
    srv_ops->close(&ds);
  }  //for (;;) //loop waiting for data to write...

  mprintf("diagsrv shutting down\n");

  diag_format_destroy(ds.format);

  srv_ops->destroy(&ds);
  mrc_domain_destroy(domain);
  free(w2);
}

void
libmrc_diag_combined_register()
{
  libmrc_diag_register_format(&ds_combined_ops);
}
