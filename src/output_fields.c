#include "psc.h"
#include "output_fields.h"
#include "util/params.h"
#include "util/profile.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

static const char *x_fldname[NR_EXTRA_FIELDS] = {
  [X_JXEX] = "jxex", [X_JYEY] = "jyey", [X_JZEZ] = "jzez",
  [X_POYX] = "poyx", [X_POYY] = "poyy", [X_POYZ] = "poyz",
  [X_E2X ] = "e2x" , [X_E2Y]  = "e2y" , [X_E2Z]  = "e2z",
  [X_B2X ] = "b2x" , [X_B2Y]  = "b2y" , [X_B2Z]  = "b2z",
};

static void
output_c_setup(struct psc_output_c *out)
{
  out->pfd.nr_flds = 0;

  fields_base_t *f = &out->pfd.flds[out->pfd.nr_flds++];
  fields_base_alloc(f, psc.ilg, psc.ihg, 3);
  f->name[0] = strdup("ne");
  f->name[1] = strdup("ni");
  f->name[2] = strdup("nn");

  f = &out->pfd.flds[out->pfd.nr_flds++];
  fields_base_alloc(f, psc.ilg, psc.ihg, 3);
  f->name[0] = strdup("jx");
  f->name[1] = strdup("jy");
  f->name[2] = strdup("jz");

  f = &out->pfd.flds[out->pfd.nr_flds++];
  fields_base_alloc(f, psc.ilg, psc.ihg, 3);
  f->name[0] = strdup("ex");
  f->name[1] = strdup("ey");
  f->name[2] = strdup("ez");

  f = &out->pfd.flds[out->pfd.nr_flds++];
  fields_base_alloc(f, psc.ilg, psc.ihg, 3);
  f->name[0] = strdup("hx");
  f->name[1] = strdup("hy");
  f->name[2] = strdup("hz");

  f = &out->pfd.flds[out->pfd.nr_flds++];
  fields_base_alloc(f, psc.ilg, psc.ihg, NR_EXTRA_FIELDS);
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    f->name[m] = strdup(x_fldname[m]);
  }

  out->tfd.nr_flds = out->pfd.nr_flds;
  for (int i = 0; i < out->pfd.nr_flds; i++) {
    fields_base_alloc(&out->tfd.flds[i], psc.ilg, psc.ihg, out->pfd.flds[i].nr_comp);
    fields_base_zero_all(&out->tfd.flds[i]);
    for (int m = 0; m < out->pfd.flds[i].nr_comp; m++) {
      out->tfd.flds[i].name[m] = strdup(out->pfd.flds[i].name[m]);
    }
  }
  out->naccum = 0;
}

#define foreach_3d(ix, iy, iz)						\
  int dx __unused = (psc.domain.ihi[0] - psc.domain.ilo[0] == 1) ? 0 : 1;	\
  int dy __unused = (psc.domain.ihi[1] - psc.domain.ilo[1] == 1) ? 0 : 1;	\
  int dz __unused = (psc.domain.ihi[2] - psc.domain.ilo[2] == 1) ? 0 : 1;	\
									\
  for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {			\
    for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {			\
      for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++)			\

#define foreach_3d_end }}

#define JX_CC(ix,iy,iz) (.5f * (F3_BASE(JXI,ix,iy,iz) + F3_BASE(JXI,ix-dx,iy,iz)))
#define JY_CC(ix,iy,iz) (.5f * (F3_BASE(JYI,ix,iy,iz) + F3_BASE(JYI,ix,iy-dy,iz)))
#define JZ_CC(ix,iy,iz) (.5f * (F3_BASE(JZI,ix,iy,iz) + F3_BASE(JZI,ix,iy,iz-dz)))

static void
calc_j(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = JX_CC(ix,iy,iz);
    XF3_BASE(f, 1, ix,iy,iz) = JY_CC(ix,iy,iz);
    XF3_BASE(f, 2, ix,iy,iz) = JZ_CC(ix,iy,iz);
  } foreach_3d_end;
}

#define EX_CC(ix,iy,iz) (.5f * (F3_BASE(EX,ix,iy,iz) + F3_BASE(EX,ix-dx,iy,iz)))
#define EY_CC(ix,iy,iz) (.5f * (F3_BASE(EY,ix,iy,iz) + F3_BASE(EY,ix,iy-dy,iz)))
#define EZ_CC(ix,iy,iz) (.5f * (F3_BASE(EZ,ix,iy,iz) + F3_BASE(EZ,ix,iy,iz-dz)))

static void
calc_E(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = EX_CC(ix,iy,iz);
    XF3_BASE(f, 1, ix,iy,iz) = EY_CC(ix,iy,iz);
    XF3_BASE(f, 2, ix,iy,iz) = EZ_CC(ix,iy,iz);
  } foreach_3d_end;
}

#define HX_CC(ix,iy,iz) (.25f*(F3_BASE(HX,ix,iy,iz   ) + F3_BASE(HX,ix,iy-dy,iz   ) + \
			       F3_BASE(HX,ix,iy,iz-dz) + F3_BASE(HX,ix,iy-dy,iz-dz)))
#define HY_CC(ix,iy,iz) (.25f*(F3_BASE(HY,ix,iy,iz   ) + F3_BASE(HY,ix-dx,iy,iz   ) + \
			       F3_BASE(HY,ix,iy,iz-dz) + F3_BASE(HY,ix-dx,iy,iz-dz)))
#define HZ_CC(ix,iy,iz) (.25f*(F3_BASE(HZ,ix,iy   ,iz) + F3_BASE(HZ,ix-dx,iy   ,iz) + \
			       F3_BASE(HZ,ix,iy-dy,iz) + F3_BASE(HZ,ix-dx,iy-dy,iz)))

static void
calc_H(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = HX_CC(ix,iy,iz);
    XF3_BASE(f, 1, ix,iy,iz) = HY_CC(ix,iy,iz);
    XF3_BASE(f, 2, ix,iy,iz) = HZ_CC(ix,iy,iz);
  } foreach_3d_end;
}

static void
output_calculate_pfields(struct psc_output_c *out)
{
  fields_base_t *p = &out->pfd.flds[out->pfd.nr_flds - 1];

  foreach_3d(ix, iy, iz) {
    XF3_BASE(p, X_JXEX, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
    XF3_BASE(p, X_JYEY, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
    XF3_BASE(p, X_JZEZ, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);

    XF3_BASE(p, X_POYX, ix,iy,iz) = (EY_CC(ix,iy,iz) * HZ_CC(ix,iy,iz) - 
				     EZ_CC(ix,iy,iz) * HY_CC(ix,iy,iz));
    XF3_BASE(p, X_POYY, ix,iy,iz) = (EZ_CC(ix,iy,iz) * HX_CC(ix,iy,iz) -
				     EX_CC(ix,iy,iz) * HZ_CC(ix,iy,iz));
    XF3_BASE(p, X_POYZ, ix,iy,iz) = (EX_CC(ix,iy,iz) * HY_CC(ix,iy,iz) -
				     EY_CC(ix,iy,iz) * HX_CC(ix,iy,iz));
    
    XF3_BASE(p, X_E2X, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
    XF3_BASE(p, X_E2Y, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
    XF3_BASE(p, X_E2Z, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
    
    XF3_BASE(p, X_B2X, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
    XF3_BASE(p, X_B2Y, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
    XF3_BASE(p, X_B2Z, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
  } foreach_3d_end;
}

static struct psc_output_format_ops *psc_output_format_ops_list[] = {
  &psc_output_format_ops_binary,
#ifdef HAVE_LIBHDF5
  &psc_output_format_ops_hdf5,
  &psc_output_format_ops_xdmf,
#endif
  &psc_output_format_ops_vtk,
  &psc_output_format_ops_vtk_points,
  &psc_output_format_ops_vtk_cells,
  NULL,
};

static struct psc_output_format_ops *
find_output_format_ops(const char *ops_name)
{
  for (int i = 0; psc_output_format_ops_list[i]; i++) {
    if (strcasecmp(psc_output_format_ops_list[i]->name, ops_name) == 0)
      return psc_output_format_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_output_format_ops '%s' not available.\n", ops_name);
  abort();
}

#define VAR(x) (void *)offsetof(struct psc_output_c, x)

static struct param psc_output_c_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")      },
  { "output_format"      , VAR(output_format)        , PARAM_STRING("binary") },
  { "output_combine"     , VAR(output_combine)       , PARAM_BOOL(0)        },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)        },
  { "pfield_first"       , VAR(pfield_first)         , PARAM_INT(0)         },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)        },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)        },
  { "tfield_first"       , VAR(tfield_first)         , PARAM_INT(0)         },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)        },

  { "output_write_jxex"  , VAR(dowrite_fd[X_JXEX])   , PARAM_BOOL(0)        },
  { "output_write_jyey"  , VAR(dowrite_fd[X_JYEY])   , PARAM_BOOL(0)        },
  { "output_write_jzez"  , VAR(dowrite_fd[X_JZEZ])   , PARAM_BOOL(0)        },
  { "output_write_poyx"  , VAR(dowrite_fd[X_POYX])   , PARAM_BOOL(0)        },
  { "output_write_poyy"  , VAR(dowrite_fd[X_POYY])   , PARAM_BOOL(0)        },
  { "output_write_poyz"  , VAR(dowrite_fd[X_POYZ])   , PARAM_BOOL(0)        },
  { "output_write_e2x"   , VAR(dowrite_fd[X_E2X])    , PARAM_BOOL(0)        },
  { "output_write_e2y"   , VAR(dowrite_fd[X_E2Y])    , PARAM_BOOL(0)        },
  { "output_write_e2z"   , VAR(dowrite_fd[X_E2Z])    , PARAM_BOOL(0)        },
  { "output_write_b2x"   , VAR(dowrite_fd[X_B2X])    , PARAM_BOOL(0)        },
  { "output_write_b2y"   , VAR(dowrite_fd[X_B2Y])    , PARAM_BOOL(0)        },
  { "output_write_b2z"   , VAR(dowrite_fd[X_B2Z])    , PARAM_BOOL(0)        },
  {},
};

#undef VAR

static struct psc_output_c psc_output_c;

// ----------------------------------------------------------------------
// output_c_create

static void output_c_create(void)
{ 
  struct psc_output_c *out = &psc_output_c;
  params_parse_cmdline(out, psc_output_c_descr, "PSC output C", MPI_COMM_WORLD);
  params_print(out, psc_output_c_descr, "PSC output C", MPI_COMM_WORLD);

  out->pfield_next = out->pfield_first;
  out->tfield_next = out->tfield_first;

  out->format_ops = find_output_format_ops(out->output_format);
  if (out->format_ops->create) {
    out->format_ops->create();
  }
};

// ----------------------------------------------------------------------
// make_fields_list

static void
make_fields_list(struct psc_fields_list *list, struct psc_fields_list *list_in,
		 bool *dowrite_fd)
{
  list->nr_flds = 0;
  for (int i = 0; i < list_in->nr_flds; i++) {
    fields_base_t *f = &list_in->flds[i];
    for (int m = 0; m < f->nr_comp; m++) {
      if (i == list_in->nr_flds - 1 && !dowrite_fd[m])
	continue;
      
      fields_base_t *fld = &list->flds[list->nr_flds++];
      fields_base_alloc_with_array(fld, psc.ilg, psc.ihg, 1,
				   &XF3_BASE(f,m, psc.ilg[0], psc.ilg[1], psc.ilg[2]));
      fld->name[0] = strdup(f->name[m]);
    }
  }
}

static void
free_fields_list(struct psc_fields_list *list)
{
  for (int m = 0; m < list->nr_flds; m++) {
    fields_base_free(&list->flds[m]);
  }
}

// ----------------------------------------------------------------------
// copy_to_global helper

static void
copy_to_global(fields_base_real_t *fld, fields_base_real_t *buf,
	       int *ilo, int *ihi, int *ilg, int *img)
{
  int *glo = psc.domain.ilo, *ghi = psc.domain.ihi;
  int my = ghi[1] - glo[1];
  int mx = ghi[0] - glo[0];

  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	fld[((iz - glo[2]) * my + iy - glo[1]) * mx + ix - glo[0]] =
	  buf[((iz - ilg[2]) * img[1] + iy - ilg[1]) * img[0] + ix - ilg[0]];
      }
    }
  }
}

// ----------------------------------------------------------------------
// write_fields_combine

char *
psc_output_c_filename(struct psc_output_c *out, const char *pfx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
  char *filename = malloc(strlen(out->data_dir) + 30);
  if (out->output_combine) {
    sprintf(filename, "%s/%s_%07d%s", out->data_dir, pfx, psc.timestep,
	    out->format_ops->ext);
  } else {
    sprintf(filename, "%s/%s_%06d_%07d%s", out->data_dir, pfx, rank, psc.timestep,
	    out->format_ops->ext);
  }
  if (rank == 0) {
    printf("[%d] write_fields: %s\n", rank, filename);
  }
  return filename;
}

// ----------------------------------------------------------------------
// write_fields_combine

static void
write_fields_combine(struct psc_output_c *out,
		     struct psc_fields_list *list, const char *prefix)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  void *ctx;
  if (rank == 0) {
    struct psc_fields_list list_combined;
    list_combined.nr_flds = list->nr_flds;

    for (int m = 0; m < list->nr_flds; m++) {
      fields_base_t *fld = &list_combined.flds[m];
      // We're creating a fake list here, with the data arrays being left
      // out. We don't want to have them all in memory at the same time,
      // and at this point, the list_combined is only used to convey info
      // about the data to come.
      fields_base_alloc_with_array(fld, psc.domain.ilo, psc.domain.ihi, 1, NULL);
      fld->name[0] = strdup(list->flds[m].name[0]);
    }
    out->format_ops->open(out, &list_combined, prefix, &ctx);
    for (int m = 0; m < list->nr_flds; m++) {
      fields_base_free(&list_combined.flds[m]);
    }
  }

  /* printf("glo %d %d %d ghi %d %d %d\n", glo[0], glo[1], glo[2], */
  /* 	     ghi[0], ghi[1], ghi[2]); */

  for (int m = 0; m < list->nr_flds; m++) {
    int s_ilo[3], s_ihi[3], s_ilg[3], s_img[3];
    fields_base_real_t *s_data = list->flds[m].flds;

    for (int d = 0; d < 3; d++) {
      s_ilo[d] = psc.ilo[d];
      s_ihi[d] = psc.ihi[d];
      s_ilg[d] = psc.ilg[d];
      s_img[d] = psc.img[d];
    }
    
    if (rank != 0) {
      MPI_Send(s_ilo, 3, MPI_INT, 0, 100, MPI_COMM_WORLD);
      MPI_Send(s_ihi, 3, MPI_INT, 0, 101, MPI_COMM_WORLD);
      MPI_Send(s_ilg, 3, MPI_INT, 0, 102, MPI_COMM_WORLD);
      MPI_Send(s_img, 3, MPI_INT, 0, 103, MPI_COMM_WORLD);
      unsigned int sz = fields_base_size(&list->flds[m]);
      MPI_Send(s_data, sz, MPI_FIELDS_BASE_REAL, 0, 104, MPI_COMM_WORLD);
    } else { // rank == 0
      fields_base_t fld;
      fields_base_alloc(&fld, psc.domain.ilo, psc.domain.ihi, 1);
      fld.name[0] = strdup(list->flds[m].name[0]);

      for (int n = 0; n < size; n++) {
	int ilo[3], ihi[3], ilg[3], img[3];
	fields_base_real_t *buf;
	
	if (n == 0) {
	  for (int d = 0; d < 3; d++) {
	    ilo[d] = s_ilo[d];
	    ihi[d] = s_ihi[d];
	    ilg[d] = s_ilg[d];
	    img[d] = s_img[d];
	  }
	  buf = s_data;
	} else {
	  MPI_Recv(ilo, 3, MPI_INT, n, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(ihi, 3, MPI_INT, n, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(ilg, 3, MPI_INT, n, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(img, 3, MPI_INT, n, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  int ntot = img[0] * img[1] * img[2];
	  buf = calloc(ntot, sizeof(*buf));
	  MPI_Recv(buf, ntot, MPI_FIELDS_BASE_REAL, n, 104, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	}
	/* printf("[%d] ilo %d %d %d ihi %d %d %d\n", rank, ilo[0], ilo[1], ilo[2], */
	/*        ihi[0], ihi[1], ihi[2]); */
	copy_to_global(fld.flds, buf, ilo, ihi, ilg, img);
	if (n != 0) {
	  free(buf);
	}
      }
      out->format_ops->write_field(ctx, &fld);
      fields_base_free(&fld);
    }
  }

  if (rank == 0) {
    out->format_ops->close(ctx);
  }

}

// ----------------------------------------------------------------------
// write_fields

static void
write_fields(struct psc_output_c *out, struct psc_fields_list *list,
	     const char *prefix)
{
  if (out->output_combine) {
    return write_fields_combine(out, list, prefix);
  }

  void *ctx;
  out->format_ops->open(out, list, prefix, &ctx);

  for (int m = 0; m < list->nr_flds; m++) {
    out->format_ops->write_field(ctx, &list->flds[m]);
  }
  
  out->format_ops->close(ctx);
}

// ----------------------------------------------------------------------
// output_c_field

static void
output_c_field()
{
  struct psc_output_c *out = &psc_output_c;

  static bool first_time = true;
  if (first_time) {
    output_c_setup(out);
    first_time = false;
  }

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  psc_calc_densities(&out->pfd.flds[0], 0); // FIXME
  calc_j(&out->pfd.flds[1]);
  calc_E(&out->pfd.flds[2]);
  calc_H(&out->pfd.flds[3]);
  output_calculate_pfields(out);

  if (out->dowrite_pfield) {
    if (psc.timestep >= out->pfield_next) {
       out->pfield_next += out->pfield_step;
       struct psc_fields_list flds_list;
       make_fields_list(&flds_list, &out->pfd, out->dowrite_fd);
       write_fields(out, &flds_list, "pfd");
       free_fields_list(&flds_list);
    }
  }

  if (out->dowrite_tfield) {
    for (int m = 0; m < out->tfd.nr_flds; m++) {
      fields_base_axpy_all(&out->tfd.flds[m], 1., &out->pfd.flds[m]); // tfd += pfd
    }
    out->naccum++;
    if (psc.timestep >= out->tfield_next) {
      out->tfield_next += out->tfield_step;

      // convert accumulated values to correct temporal mean
      for (int m = 0; m < out->tfd.nr_flds; m++) {
	fields_base_scale_all(&out->tfd.flds[m], 1. / out->naccum);
      }

      struct psc_fields_list flds_list;
      make_fields_list(&flds_list, &out->tfd, out->dowrite_fd);
      write_fields(out, &flds_list, "tfd");
      free_fields_list(&flds_list);
      for (int m = 0; m < out->tfd.nr_flds; m++) {
	fields_base_zero_all(&out->tfd.flds[m]);
      }
      out->naccum = 0;
    }
  }
  
  prof_stop(pr);
}

// ======================================================================
// psc_output_ops_c

struct psc_output_ops psc_output_ops_c = {
  .name           = "c",
  .create         = output_c_create,
  .out_field      = output_c_field,
};

