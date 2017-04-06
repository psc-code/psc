
#include "ggcm_mhd_diag_private.h"
#include "ggcm_mhd_diag_item_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd.h"

#include <mrc_io.h>
#include <mrc_fld_as_float.h>

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// ======================================================================
// ggcm_mhd_diag_c

struct mrc_io_entry {
  struct mrc_io *io;
  int outtype;
  float sheet;
  list_t entry;
};

static struct mrc_io *
create_mrc_io(struct ggcm_mhd_diag *diag, int outtype, float sheet)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);
  struct mrc_io *io;

  if (sub->rank_diagsrv > 0) {
    io = ggcm_diag_lib_create_mrc_io(ggcm_mhd_diag_comm(diag), sub->run, "combined",
				     outtype, sheet, sub->rank_diagsrv);
  } else {
    const char *outputmode = "xdmf_collective";
    mrc_params_get_option_string("outputmode", &outputmode);

    io = ggcm_diag_lib_create_mrc_io(ggcm_mhd_diag_comm(diag), sub->run, outputmode,
				     outtype, sheet, -1);
  }

  return io;
}

static struct mrc_io_entry *
find_mrc_io(struct ggcm_mhd_diag *diag, int outtype, float sheet)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);
  
  struct mrc_io_entry *p;
  list_for_each_entry(p, &sub->mrc_io_list, entry) {
    if (p->outtype == outtype && p->sheet == sheet) {
      return p;
    }
  }
  return NULL;
}

static struct mrc_io_entry *
get_mrc_io(struct ggcm_mhd_diag *diag, int outtype, float sheet)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);

  struct mrc_io_entry *p = find_mrc_io(diag, outtype, sheet);
  if (p)
    return p;

  p = calloc(1, sizeof(*p));
  p->io = create_mrc_io(diag, outtype, sheet);
  p->outtype = outtype;
  p->sheet = sheet;
  list_add_tail(&p->entry, &sub->mrc_io_list);
  return p;
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_create

static void
ggcm_mhd_diag_c_create(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);

  INIT_LIST_HEAD(&sub->mrc_io_list);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_setup_common
//
// shared between setup() and read(), though really, we should 
// restore from checkpoint for read(), rather than parsing stuff again

static void
ggcm_mhd_diag_c_setup_common(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);

  sub->nr_planes[0] = 
    parse_float_array(sub->outplanex, sub->planes[0], MAX_PLANES);
  sub->nr_planes[1] = 
    parse_float_array(sub->outplaney, sub->planes[1], MAX_PLANES);
  sub->nr_planes[2] = 
    parse_float_array(sub->outplanez, sub->planes[2], MAX_PLANES);

  // parse comma/colon separated list of fields
  char *s_orig = strdup(sub->fields), *p, *s = s_orig;
  while ((p = strsep(&s, ":, "))) {
    struct ggcm_mhd_diag_item *item =
      ggcm_mhd_diag_item_create(ggcm_mhd_diag_comm(diag));
    ggcm_mhd_diag_item_set_type(item, p);
    ggcm_mhd_diag_item_set_param_obj(item, "diag", diag);
    ggcm_mhd_diag_add_child(diag, (struct mrc_obj *) item);
  }
  free(s_orig);

}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_setup

static void
ggcm_mhd_diag_c_setup(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);

  if (sub->find_diagsrv) {
    sub->rank_diagsrv = sub->find_diagsrv(diag);
  } else {
    sub->rank_diagsrv = -1;
  }

  ggcm_mhd_diag_c_setup_common(diag);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_read

static void
ggcm_mhd_diag_c_read(struct ggcm_mhd_diag *diag, struct mrc_io *io)
{
  ggcm_mhd_diag_c_create(diag);
  ggcm_mhd_diag_c_setup_common(diag);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_write_one_fld

void
ggcm_mhd_diag_c_write_one_fld(struct mrc_io *io, struct mrc_fld *_fld,
			      int outtype, float plane)
{
  // FIXME, this works around mrc_fld_write() not handling AOS layout fields
  struct mrc_fld *fld = mrc_fld_get_as(_fld, FLD_TYPE);

  switch (outtype) {
  case DIAG_TYPE_3D:
    mrc_fld_write(fld, io);
    break;
  case DIAG_TYPE_2D_X:
  case DIAG_TYPE_2D_Y:
  case DIAG_TYPE_2D_Z:
    mrc_io_write_field_slice(io, 1., fld, outtype, plane);
    break;
  default:
    assert(0);
  }

  mrc_fld_put_as(fld, _fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_write_one_field

void
ggcm_mhd_diag_c_write_one_field(struct mrc_io *io, struct mrc_fld *_f, int m,
				const char *name, float scale, int outtype,
				float plane)
{
  struct mrc_fld *f = mrc_fld_get_as(_f, FLD_TYPE);

  assert(m < f->_nr_comps);

  int bnd = f->_nr_ghosts;
  struct mrc_fld *fld = mrc_domain_fld_create(f->_domain, bnd, name);
  mrc_fld_set_type(fld, FLD_TYPE);
  char s[strlen(name) + 10];
  sprintf(s, "mrc_fld_%s", name);
  mrc_fld_set_name(fld, s);
  mrc_fld_setup(fld);

  // ghosts may not be set
  for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
    mrc_fld_foreach(fld, ix,iy,iz, bnd, bnd) {
      M3(fld,0, ix,iy,iz, p) = scale * M3(f,m, ix,iy,iz, p);
    } mrc_fld_foreach_end;
  }

  ggcm_mhd_diag_c_write_one_fld(io, fld, outtype, plane);

  mrc_fld_destroy(fld);

  // FIXME, should use _put_as, but don't want copy-back
  if (strcmp(mrc_fld_type(_f), FLD_TYPE) != 0) {
    mrc_fld_destroy(f);
  }
}

// ----------------------------------------------------------------------
// write_fields

static void
write_fields(struct ggcm_mhd_diag *diag, struct mrc_fld *fld,
	     struct mrc_io *io, int diag_type, float plane)
{
  struct ggcm_mhd_diag_item *item;
  list_for_each_entry(item, &diag->obj.children_list, obj.child_entry) {
    ggcm_mhd_diag_item_run(item, io, fld, diag_type, plane);
  }

#if 0
  if(use_rcm.and.rcm_diag) {
    diag_write_field_c(scrr,rcmrrl,#nnx,#nny,#nnz,'rcmrr',diag_type_3d,-1.,-99);
    diag_write_field_c(scpp,rcmppl,#nnx,#nny,#nnz,'rcmpp',diag_type_3d,-1.,-99);
    diag_write_field_c(1.,fbmask,#nnx,#nny,#nnz,'rcmmask',diag_type_3d,-1.,-99);
  }
#endif
}

// ----------------------------------------------------------------------
// diagc3

static void
diagc3(struct ggcm_mhd_diag *diag, struct mrc_fld *fld, int itdia,
       char *time_str)
{
  struct ggcm_mhd *mhd = diag->mhd;

  struct mrc_io *io = get_mrc_io(diag, DIAG_TYPE_3D, -1)->io;
  mrc_io_open(io, "w", itdia, mhd->time_code * mhd->tnorm);
  ggcm_diag_lib_write_openggcm_attrs(io, time_str);
  write_fields(diag, fld, io, DIAG_TYPE_3D, -1);
  mrc_io_close(io);
}

// ----------------------------------------------------------------------
// diagcxyz

static void
diagcxyz(struct ggcm_mhd_diag *diag, struct mrc_fld *fld, int itdia,
	 char *time_str, int d)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);
  struct ggcm_mhd *mhd = diag->mhd;

  int diag_type = DIAG_TYPE_2D_X + d;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  const double *lo = mrc_crds_lo(crds), *hi = mrc_crds_hi(crds);

  for (int i = 0; i < sub->nr_planes[d]; i++) {
    float plane = sub->planes[d][i];
    if (plane < lo[d] || plane > hi[d])
      continue;

    struct mrc_io *io = get_mrc_io(diag, diag_type, plane)->io;
    mrc_io_open(io, "w", itdia, mhd->time_code * mhd->tnorm);
    ggcm_diag_lib_write_openggcm_attrs(io, time_str);
    write_fields(diag, fld, io, diag_type, plane);
    mrc_io_close(io);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_run_now

static void
ggcm_mhd_diag_c_run_now(struct ggcm_mhd_diag *diag, struct mrc_fld *fld,
			int diag_type, int itdia)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);
  struct ggcm_mhd *mhd = diag->mhd;

  char time_str[80] = "TIME";
  if (sub->make_time_string) {
    sub->make_time_string(time_str, mhd->time_code * mhd->tnorm, mhd->dacttime);
  }

  switch (diag_type) {
  case DIAG_TYPE_2D_X:
    mpi_printf(ggcm_mhd_diag_comm(diag),
	     "ggcm_mhd_diag: writing 2d output (step %d)\n", itdia);
    /* fall through */
  case DIAG_TYPE_2D_Y:
  case DIAG_TYPE_2D_Z:
    diagcxyz(diag, fld, itdia, time_str, diag_type - DIAG_TYPE_2D_X);
    break;

  case DIAG_TYPE_3D:
    mpi_printf(ggcm_mhd_diag_comm(diag),
	     "ggcm_mhd_diag: writing 3d output (step %d)\n", itdia);
    diagc3(diag, fld, itdia, time_str);
    break;

  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_run

static void
ggcm_mhd_diag_c_run(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);
  struct ggcm_mhd *mhd = diag->mhd;
  bool output3d, output2d;
  int itdia3d, itdia2d;

  // FIXME, ggcm_mhd_diag_run() is an entirely OpenGGCM-specific interface only
  // which really should be more generic...
  if (sub->run_hack) {
    sub->run_hack(diag, &output3d, &output2d, &itdia3d, &itdia2d);
  } else {
    output3d = true;
    output2d = false;
    itdia3d = mhd->istep;
    itdia2d = mhd->istep;
  }

  if (!(output2d || output3d))
    return;

  ggcm_mhd_fill_ghosts(mhd, mhd->fld, mhd->time_code);

  if (output3d) {
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_3D, itdia3d);
  }

  if (output2d) {
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_2D_X + 0, itdia2d);
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_2D_X + 1, itdia2d);
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_2D_X + 2, itdia2d);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_shutdown

static void
ggcm_mhd_diag_c_shutdown(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);

  while(!list_empty(&sub->mrc_io_list)) {
    struct mrc_io_entry *p =
      list_entry(sub->mrc_io_list.next, struct mrc_io_entry, entry);
    mrc_io_destroy(p->io);
    list_del(&p->entry);
    free(p);
  }
}

// ----------------------------------------------------------------------
// diagsrv_one_c

static void
diagsrv_one_c(struct mrc_mod *mod, void *arg)
{
  int rc;

  const char *outputmode, *diagsc_srv;
  rc = mrc_params_get_option_string("outputmode", &outputmode);
  assert(rc == 0);
  rc = mrc_params_get_option_string("diagsc_srv", &diagsc_srv);
  assert(rc == 0);
  int n_mhd_procs = (unsigned long) arg;

  mrc_io_server(outputmode, diagsc_srv, n_mhd_procs);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_mod_register

static void
ggcm_mhd_diag_c_mod_register(struct ggcm_mhd_diag *mhd_diag, struct mrc_mod *mod)
{
  const char *outputmode = NULL;
  mrc_params_get_option_string("outputmode", &outputmode);

  if (strcmp(outputmode, "xdmf") != 0 &&
      strcmp(outputmode, "xdmf_collective") != 0) {
    int n_mhd_procs = mrc_mod_get_nr_procs(mod, "MHD");
    mrc_mod_register(mod, "DIAGSC", 1, diagsrv_one_c, (void *)(unsigned long) n_mhd_procs);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag c subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_diag_c, x)
static struct param ggcm_mhd_diag_c_descr[] = {
  { "run"             , VAR(run)             , PARAM_STRING("run")        },
  { "fields"          , VAR(fields)          , PARAM_STRING("rr:pp:v:b")  },
  { "outplanex"       , VAR(outplanex)       , PARAM_STRING("")           },
  { "outplaney"       , VAR(outplaney)       , PARAM_STRING("")           },
  { "outplanez"       , VAR(outplanez)       , PARAM_STRING("")           },

  { "rank_diagsrv"    , VAR(rank_diagsrv)    , MRC_VAR_INT                },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_diag subclass "c"

struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops = {
  .name         = "c",
  .size         = sizeof(struct ggcm_mhd_diag_c),
  .param_descr  = ggcm_mhd_diag_c_descr,
  .create       = ggcm_mhd_diag_c_create,
  .setup        = ggcm_mhd_diag_c_setup,
  .read         = ggcm_mhd_diag_c_read,
  .run          = ggcm_mhd_diag_c_run,
  .run_now      = ggcm_mhd_diag_c_run_now,
  .shutdown     = ggcm_mhd_diag_c_shutdown,
  .mod_register = ggcm_mhd_diag_c_mod_register,
};
