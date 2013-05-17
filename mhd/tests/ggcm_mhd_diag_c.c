
#include "ggcm_mhd_diag_private.h"
#include "ggcm_mhd_diag_item_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd.h"

#include <mrc_io.h>
#include <mrc_mod.h>

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

// ======================================================================
// ggcm_mhd_diag_c

#define MAX_PLANES (10)

struct ggcm_mhd_diag_c {
  // parameters
  char *fields;
  char *outplanex;
  char *outplaney;
  char *outplanez;
  float planes[3][MAX_PLANES];
  float nr_planes[3];

  // state
  int rank_diagsrv;
  list_t mrc_io_list;
};

#define ggcm_mhd_diag_c(diag) mrc_to_subobj(diag, struct ggcm_mhd_diag_c)

// ----------------------------------------------------------------------

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
    io = ggcm_diag_lib_create_mrc_io(ggcm_mhd_diag_comm(diag), "combined",
				     outtype, sheet, sub->rank_diagsrv);
  } else {
    const char *outputmode = "xdmf_collective";
    mrc_params_get_option_string("outputmode", &outputmode);

    io = ggcm_diag_lib_create_mrc_io(ggcm_mhd_diag_comm(diag), outputmode,
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
// ggcm_mhd_diag_c_setup

static void
ggcm_mhd_diag_c_setup(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);

  int size;
  MPI_Comm_size(ggcm_mhd_diag_comm(diag), &size);
  //  if (size == 1 || !_ggcm) {
  if (true) {
    sub->rank_diagsrv = -1;
  } else {
    //    sub->rank_diagsrv = mrc_mod_get_first_node(_ggcm->mod, "DIAGSC");
  }

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
// ggcm_mhd_diag_c_read

static void
ggcm_mhd_diag_c_read(struct ggcm_mhd_diag *diag, struct mrc_io *io)
{
  ggcm_mhd_diag_c_create(diag);
  ggcm_mhd_diag_c_setup(diag);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_write_one_f3

void
ggcm_mhd_diag_c_write_one_f3(struct mrc_io *io, struct mrc_f3 *fld,
			     int outtype, float plane)
{
  switch (outtype) {
  case DIAG_TYPE_3D: {
    mrc_f3_write(fld, io);
    break;
  }
  case DIAG_TYPE_2D_X:
  case DIAG_TYPE_2D_Y:
  case DIAG_TYPE_2D_Z: {
    mrc_io_write_field_slice(io, 1., fld, outtype, plane);
    break;
  }
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_write_one_field

void
ggcm_mhd_diag_c_write_one_field(struct mrc_io *io, struct mrc_f3 *f, int m,
				const char *name, float scale, int outtype,
				float plane)
{
  struct mrc_f3 *fld = mrc_domain_f3_create(f->_domain, SW_2, name);
  char s[strlen(name) + 10];
  sprintf(s, "mrc_f3_%s", name);
  mrc_f3_set_name(fld, s);
  mrc_f3_setup(fld);
  mrc_f3_foreach(fld, ix,iy,iz, 2, 2) {
    MRC_F3(fld,0, ix,iy,iz) = scale * MRC_F3(f,m, ix,iy,iz);
  } mrc_f3_foreach_end;

  ggcm_mhd_diag_c_write_one_f3(io, fld, outtype, plane);

  mrc_f3_destroy(fld);
}

// ----------------------------------------------------------------------
// write_fields

static void
write_fields(struct ggcm_mhd_diag *diag, struct mrc_io *io, int diag_type, float plane)
{
  struct ggcm_mhd *mhd = diag->mhd;
  struct ggcm_mhd_flds *flds = ggcm_mhd_flds_get_as(mhd->flds_base, "c");
  struct mrc_f3 *f = ggcm_mhd_flds_get_mrc_f3(flds);

  struct ggcm_mhd_diag_item *item;
  list_for_each_entry(item, &diag->obj.children_list, obj.child_entry) {
    ggcm_mhd_diag_item_run(item, io, f, diag_type, plane);
  }

#if 0
  if(use_rcm.and.rcm_diag) {
    diag_write_field_c(scrr,rcmrrl,#nnx,#nny,#nnz,'rcmrr',diag_type_3d,-1.,-99);
    diag_write_field_c(scpp,rcmppl,#nnx,#nny,#nnz,'rcmpp',diag_type_3d,-1.,-99);
    diag_write_field_c(1.,fbmask,#nnx,#nny,#nnz,'rcmmask',diag_type_3d,-1.,-99);
  }
#endif

  ggcm_mhd_flds_put_as(flds, mhd->flds_base);
}

// ----------------------------------------------------------------------
// diagc3

static void
diagc3(struct ggcm_mhd_diag *diag, int itdia, char *time_str)
{
  struct ggcm_mhd *mhd = diag->mhd;

  struct mrc_io *io = get_mrc_io(diag, DIAG_TYPE_3D, -1)->io;
  mrc_io_open(io, "w", itdia, mhd->time);
  ggcm_diag_lib_write_openggcm_attrs(io, time_str);
  write_fields(diag, io, DIAG_TYPE_3D, -1);
  mrc_io_close(io);
}

// ----------------------------------------------------------------------
// diagcxyz

static void
diagcxyz(struct ggcm_mhd_diag *diag, int itdia, char *time_str, int d)
{
  struct ggcm_mhd_diag_c *sub = ggcm_mhd_diag_c(diag);
  struct ggcm_mhd *mhd = diag->mhd;

  int diag_type = DIAG_TYPE_2D_X + d;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  float xl[3], xh[3];
  mrc_crds_get_param_float3(crds, "l", xl);
  mrc_crds_get_param_float3(crds, "h", xh);

  for (int i = 0; i < sub->nr_planes[d]; i++) {
    float plane = sub->planes[d][i];
    if (plane < xl[d] || plane > xh[d])
      continue;

    struct mrc_io *io = get_mrc_io(diag, diag_type, plane)->io;
    mrc_io_open(io, "w", itdia, mhd->time);
    ggcm_diag_lib_write_openggcm_attrs(io, time_str);
    write_fields(diag, io, diag_type, plane);
    mrc_io_close(io);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_c_run

static void
ggcm_mhd_diag_c_run(struct ggcm_mhd_diag *diag)
{
  struct ggcm_mhd *mhd = diag->mhd;
  int itdia3d, itdia2d;

  bool output3d = true;
  bool output2d = false;
  itdia3d = mhd->istep;
  itdia2d = mhd->istep;

  ggcm_mhd_fill_ghosts(mhd, _RR1, mhd->time);

  char time_str[80] = "TIME";
  //  ggcm_diag_lib_make_time_string(time_str, mhd->time, mhd->dacttime);

  if (output3d) {
    mpi_printf(ggcm_mhd_diag_comm(diag),
	       "ggcm_mhd_diag: writing 3d output (step %d)\n", itdia3d);
    diagc3(diag, itdia3d, time_str);
  }

  if (output2d) {
    mpi_printf(ggcm_mhd_diag_comm(diag),
	       "ggcm_mhd_diag: writing 2d output (step %d)\n", itdia2d);
    diagcxyz(diag, itdia2d, time_str, 0);
    diagcxyz(diag, itdia2d, time_str, 1);
    diagcxyz(diag, itdia2d, time_str, 2);
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
// ggcm_mhd_diag c subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_diag_c, x)
static struct param ggcm_mhd_diag_c_descr[] = {
  { "fields"          , VAR(fields)          , PARAM_STRING("rr:pp:v:b")  },
  { "outplanex"       , VAR(outplanex)       , PARAM_STRING("")           },
  { "outplaney"       , VAR(outplaney)       , PARAM_STRING("")           },
  { "outplanez"       , VAR(outplanez)       , PARAM_STRING("")           },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_diag subclass "c"

struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops = {
  .name        = "c",
  .size        = sizeof(struct ggcm_mhd_diag_c),
  .param_descr = ggcm_mhd_diag_c_descr,
  .create      = ggcm_mhd_diag_c_create,
  .setup       = ggcm_mhd_diag_c_setup,
  .read        = ggcm_mhd_diag_c_read,
  .run         = ggcm_mhd_diag_c_run,
  .shutdown    = ggcm_mhd_diag_c_shutdown,
};

// ======================================================================

// ----------------------------------------------------------------------
// parse_float_array
//
// parse colon-separated array of floats, up to n elements
// returns the actual number of elements found in the string

int
parse_float_array(const char *str, float *arr, int n)
{
  char *s1, *s = strdup(str);
  int i;
  for (i = 0; (s1 = strsep(&s, ":")) && i < n; i++) {
    arr[i] = atof(s1);
  }
  free(s);
  return i;
}

