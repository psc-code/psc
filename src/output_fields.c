#include "psc.h"
#include "output_fields.h"
#include "util/params.h"
#include "util/profile.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

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
calc_jdote(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
    XF3_BASE(f, 1, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
    XF3_BASE(f, 2, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);
  } foreach_3d_end;
}

static void
calc_poyn(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = (EY_CC(ix,iy,iz) * HZ_CC(ix,iy,iz) - 
				EZ_CC(ix,iy,iz) * HY_CC(ix,iy,iz));
    XF3_BASE(f, 1, ix,iy,iz) = (EZ_CC(ix,iy,iz) * HX_CC(ix,iy,iz) -
				EX_CC(ix,iy,iz) * HZ_CC(ix,iy,iz));
    XF3_BASE(f, 2, ix,iy,iz) = (EX_CC(ix,iy,iz) * HY_CC(ix,iy,iz) -
				EY_CC(ix,iy,iz) * HX_CC(ix,iy,iz));
  } foreach_3d_end;
}

static void
calc_E2(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
    XF3_BASE(f, 1, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
    XF3_BASE(f, 2, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
  } foreach_3d_end;
}

static void
calc_H2(fields_base_t *f)
{
  foreach_3d(ix, iy, iz) {
    XF3_BASE(f, 0, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
    XF3_BASE(f, 1, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
    XF3_BASE(f, 2, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
  } foreach_3d_end;
}

struct output_field {
  char *name;
  int nr_comp;
  char *fld_names[3];
  void (*calc)(fields_base_t *f);
};

static struct output_field output_fields[] = {
  { .name = "n"    , .nr_comp = 3, .fld_names = { "ne", "ni", "nn" },
    .calc = psc_calc_densities },
  { .name = "j"    , .nr_comp = 3, .fld_names = { "jx", "jy", "jz" },
    .calc = calc_j },
  { .name = "e"    , .nr_comp = 3, .fld_names = { "ex", "ey", "ez" },
    .calc = calc_E },
  { .name = "h"    , .nr_comp = 3, .fld_names = { "hx", "hy", "hz" },
    .calc = calc_H },
  { .name = "jdote", .nr_comp = 3, .fld_names = { "jxex", "jyey", "jzez" },
    .calc = calc_jdote },
  { .name = "poyn" , .nr_comp = 3, .fld_names = { "poynx", "poyny", "poynz" },
    .calc = calc_poyn },
  { .name = "e2"   , .nr_comp = 3, .fld_names = { "ex2", "ey2", "ez2" },
    .calc = calc_E2 },
  { .name = "h2"   , .nr_comp = 3, .fld_names = { "hx2", "hy2", "hz2" },
    .calc = calc_H2 },
  {},
};

static struct output_field *
find_output_field(const char *name)
{
  for (int i = 0; output_fields[i].name; i++) {
    struct output_field *of = &output_fields[i];
    if (strcasecmp(of->name, name) == 0) {
      return of;
    }
  }
  fprintf(stderr, "ERROR: output_field '%s' unknown!\n", name);
  abort();
}

static void
output_c_setup(struct psc_output_c *out)
{
  struct psc_fields_list *pfd = &out->pfd;

  // setup pfd according to output_fields as given
  // (potentially) on the command line
  pfd->nr_flds = 0;
  // parse comma separated list of fields
  char *p, *s = strdup(out->output_fields);
  while ((p = strsep(&s, ", "))) {
    struct output_field *of = find_output_field(p);
    fields_base_t *f = &pfd->flds[pfd->nr_flds];
    out->out_flds[pfd->nr_flds] = of;
    pfd->nr_flds++;
      
    fields_base_alloc(f, psc.ilg, psc.ihg, of->nr_comp);
    for (int m = 0; m < of->nr_comp; m++) {
      f->name[m] = strdup(of->fld_names[m]);
    }
  }
  free(s);

  // create tfd to look just like pfd
  // FIXME, only if necessary
  struct psc_fields_list *tfd = &out->tfd;
  tfd->nr_flds = pfd->nr_flds;
  for (int i = 0; i < pfd->nr_flds; i++) {
    fields_base_alloc(&tfd->flds[i], psc.ilg, psc.ihg, pfd->flds[i].nr_comp);
    fields_base_zero_all(&tfd->flds[i]);
    for (int m = 0; m < pfd->flds[i].nr_comp; m++) {
      tfd->flds[i].name[m] = strdup(pfd->flds[i].name[m]);
    }
  }
  out->naccum = 0;
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
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "output_format"      , VAR(output_format)        , PARAM_STRING("binary")  },
  { "output_combine"     , VAR(output_combine)       , PARAM_BOOL(0)           },
  { "output_fields"      , VAR(output_fields)        , PARAM_STRING("n,j,e,h") },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)           },
  { "pfield_first"       , VAR(pfield_first)         , PARAM_INT(0)            },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)           },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)           },
  { "tfield_first"       , VAR(tfield_first)         , PARAM_INT(0)            },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)           },
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
make_fields_list(struct psc_fields_list *list, struct psc_fields_list *list_in)
{
  // the only thing this still does is to flatten
  // the list so that it only contains 1-component entries

  list->nr_flds = 0;
  for (int i = 0; i < list_in->nr_flds; i++) {
    fields_base_t *f = &list_in->flds[i];
    for (int m = 0; m < f->nr_comp; m++) {
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
// psc_output_c_filename

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
    out->format_ops->open(out, NULL, prefix, &ctx);
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

  if ((out->dowrite_pfield && psc.timestep >= out->pfield_next) ||
      out->dowrite_tfield) {
    struct psc_fields_list *pfd = &out->pfd;
    for (int i = 0; i < pfd->nr_flds; i++) {
      out->out_flds[i]->calc(&pfd->flds[i]);
    }
  }
  
  if (out->dowrite_pfield) {
    if (psc.timestep >= out->pfield_next) {
       out->pfield_next += out->pfield_step;
       struct psc_fields_list flds_list;
       make_fields_list(&flds_list, &out->pfd);
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
      make_fields_list(&flds_list, &out->tfd);
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

