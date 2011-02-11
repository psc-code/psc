
#include "mrc_diag_private.h"
#include "mrc_params.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

// ======================================================================
// diag_format

#define check_is_setup(format) do { assert(format->is_setup); } while (0)

static inline struct diag_format_ops *
diag_format_ops(struct diag_format *format)
{
  return (struct diag_format_ops *) format->obj.ops;
}

// ----------------------------------------------------------------------
// diag_format_setup

static void
_diag_format_setup(struct mrc_obj *obj)
{
  struct diag_format *format = to_diag_format(obj);
  if (!format->diag_info.outdir) {
    format->diag_info.outdir = ".";
  }

  MPI_Comm_rank(format->obj.comm, &format->rank);
  if (diag_format_ops(format)->setup) {
    diag_format_ops(format)->setup(&format->obj);
  }
  format->is_setup = true;
}

// ----------------------------------------------------------------------
// diag_format_open

void
diag_format_open(struct diag_format *format, float sheet, int outtype, int step, float time,
		 const char *time_str)
{
  check_is_setup(format);
  struct diag_format_ops *ops = diag_format_ops(format);
  format->time = time;
  format->time_str = strdup(time_str);
  ops->open(format, sheet, outtype, step);
}

// ----------------------------------------------------------------------
// diag_format_close

void
diag_format_close(struct diag_format *format)
{
  struct diag_format_ops *ops = diag_format_ops(format);
  ops->close(format);
  free(format->time_str);
}

// ----------------------------------------------------------------------
// diag_format_write_field

void
diag_format_write_field(struct diag_format *format, float scale, struct mrc_f3 *fld, int m)
{
  struct diag_format_ops *ops = diag_format_ops(format);
  if (!fld->name[m]) {
    char s[10];
    sprintf(s, "m%d", m);
    fld->name[m] = strdup(s);
  }
  ops->write_field(format, scale, fld, m);
}

// ----------------------------------------------------------------------
// diag_format_write_field2d

void
diag_format_write_field2d(struct diag_format *format, float scale, struct mrc_f2 *fld,
			  const char *fld_name, int outtype, float sheet)
{
  struct diag_format_ops *ops = diag_format_ops(format);
  ops->write_field2d(format, scale, fld, fld_name, outtype, sheet);
}

// ----------------------------------------------------------------------
// diag_format_write_field_slice

void
diag_format_write_field_slice(struct diag_format *format, float scale, struct mrc_f3 *fld,
			      const char *fld_name, int outtype, float sheet)
{
  struct diag_format_ops *ops = diag_format_ops(format);

  //dig out a sheet of constant x/y/z

  assert(outtype >= DIAG_TYPE_2D_X && outtype <= DIAG_TYPE_2D_Z);
  int dim = outtype - DIAG_TYPE_2D_X;

  int dims[3];
  mrc_domain_get_local_offset_dims(fld->domain, NULL, dims);

  //check for existence on local proc.
  //0,1, nnx-2, nnx-1 are the ghostpoints
  struct mrc_crds *crds = mrc_domain_get_crds(fld->domain);
  float *crd = crds->crd[dim];
  if (crd[1] < sheet && sheet <= crd[dims[dim]+1]) { 
    int ii;
    for (ii = 1; ii < dims[dim]; ii++) {
      if (sheet <= crd[ii+1])
	break;
    }

    float s2 = (sheet - crd[ii]) / (crd[ii+1] - crd[ii]);
    float s1 = 1. - s2;
    s1 *= scale;
    s2 *= scale;

    struct mrc_f2 f2;
    switch (dim) {
    case 0:
      mrc_f2_alloc(&f2, NULL, (int [2]) { dims[1], dims[2] }, 1);
      for(int iz = 0; iz < dims[2]; iz++) {
	for(int iy = 0; iy < dims[1]; iy++) {
	  MRC_F2(&f2,0, iy,iz) = (s1 * MRC_F3(fld,0, ii  ,iy+2,iz+2) +
				  s2 * MRC_F3(fld,0, ii+1,iy+2,iz+2));
	}
      }
      break;
    case 1:
      mrc_f2_alloc(&f2, NULL, (int [2]) { dims[0], dims[2] }, 1);
      for(int iz = 0; iz < dims[2]; iz++) {
	for(int ix = 0; ix < dims[0]; ix++) {
	  MRC_F2(&f2,0, ix,iz) = (s1 * MRC_F3(fld,0, ix+2,ii  ,iz+2) +
				  s2 * MRC_F3(fld,0, ix+2,ii+1,iz+2));
	}
      }
      break;
    case 2:
      mrc_f2_alloc(&f2, NULL, (int [2]) { dims[0], dims[1] }, 1);
      for(int iy = 0; iy < dims[1]; iy++) {
	for(int ix = 0; ix < dims[0]; ix++) {
	  MRC_F2(&f2,0, ix,iy) = (s1 * MRC_F3(fld,0, ix+2,iy+2,ii  ) +
				  s2 * MRC_F3(fld,0, ix+2,iy+2,ii+1));
	}
      }
      break;
    }

    f2.domain = fld->domain;
    ops->write_field2d(format, 1., &f2, fld_name, outtype, sheet);
    mrc_f2_free(&f2);
  } else {
    struct mrc_f2 f2 = {};
    f2.domain = fld->domain;
    ops->write_field2d(format, 1., &f2, fld_name, outtype, sheet);
  }
}

// ----------------------------------------------------------------------
// diag_format_write_fields

void
diag_format_write_fields(struct diag_format *format, struct mrc_f3 *f, int mm[])
{
  if (!mm) { // all components
    for (int m = 0; m < f->nr_comp; m++) {
      diag_format_write_field(format, 1.f, f, m);
    }
  } else {
    for (int i = 0; mm[i] >= 0; i++) {
      diag_format_write_field(format, 1.f, f, mm[i]);
    }
  }
}

// ======================================================================
// diag_format class

static LIST_HEAD(diag_formats);

void
libmrc_diag_register_format(struct diag_format_ops *ops)
{
  list_add_tail(&ops->list, &diag_formats);
}

// ----------------------------------------------------------------------
// diag_format_init

static void
diag_format_init()
{
  libmrc_diag_combined_register();
  libmrc_diag_ascii_register();
#ifdef HAVE_HDF5_H
  libmrc_diag_hdf5_register();
#endif
}

#define VAR(x) (void *)offsetof(struct diag_info, x)
static struct param diag_info_descr[] = {
  { "outdir"          , VAR(outdir)       , PARAM_STRING(".")      },
  { "run"             , VAR(run)          , PARAM_STRING("run")    },
  {},
};
#undef VAR

struct mrc_class mrc_class_diag_format = {
  .name         = "diag_format",
  .size         = sizeof(struct diag_format),
  .subclasses   = &diag_formats,
  .param_descr  = diag_info_descr,
  .param_offset = offsetof(struct diag_format, diag_info),
  .init         = diag_format_init,
  .setup        = _diag_format_setup,
};

