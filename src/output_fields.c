#include "psc.h"
#include "output_fields.h"
#include <mrc_params.h>
#include <mrc_profile.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define define_dxdydz(dx, dy, dz)					\
  int dx __unused = (psc.domain.gdims[0] == 1) ? 0 : 1;			\
  int dy __unused = (psc.domain.gdims[1] == 1) ? 0 : 1;			\
  int dz __unused = (psc.domain.gdims[2] == 1) ? 0 : 1

#define JX_CC(ix,iy,iz) (.5f * (F3_BASE(&psc.pf, JXI,ix,iy,iz) + F3_BASE(&psc.pf, JXI,ix-dx,iy,iz)))
#define JY_CC(ix,iy,iz) (.5f * (F3_BASE(&psc.pf, JYI,ix,iy,iz) + F3_BASE(&psc.pf, JYI,ix,iy-dy,iz)))
#define JZ_CC(ix,iy,iz) (.5f * (F3_BASE(&psc.pf, JZI,ix,iy,iz) + F3_BASE(&psc.pf, JZI,ix,iy,iz-dz)))

static void
calc_j(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = JX_CC(ix,iy,iz);
      F3_BASE(f, 1, ix,iy,iz) = JY_CC(ix,iy,iz);
      F3_BASE(f, 2, ix,iy,iz) = JZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

#define EX_CC(ix,iy,iz) (.5f * (F3_BASE(&psc.pf, EX,ix,iy,iz) + F3_BASE(&psc.pf, EX,ix-dx,iy,iz)))
#define EY_CC(ix,iy,iz) (.5f * (F3_BASE(&psc.pf, EY,ix,iy,iz) + F3_BASE(&psc.pf, EY,ix,iy-dy,iz)))
#define EZ_CC(ix,iy,iz) (.5f * (F3_BASE(&psc.pf, EZ,ix,iy,iz) + F3_BASE(&psc.pf, EZ,ix,iy,iz-dz)))

static void
calc_E(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = EX_CC(ix,iy,iz);
      F3_BASE(f, 1, ix,iy,iz) = EY_CC(ix,iy,iz);
      F3_BASE(f, 2, ix,iy,iz) = EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

#define HX_CC(ix,iy,iz) (.25f*(F3_BASE(&psc.pf, HX,ix,iy,iz   ) + F3_BASE(&psc.pf, HX,ix,iy-dy,iz   ) + \
			       F3_BASE(&psc.pf, HX,ix,iy,iz-dz) + F3_BASE(&psc.pf, HX,ix,iy-dy,iz-dz)))
#define HY_CC(ix,iy,iz) (.25f*(F3_BASE(&psc.pf, HY,ix,iy,iz   ) + F3_BASE(&psc.pf, HY,ix-dx,iy,iz   ) + \
			       F3_BASE(&psc.pf, HY,ix,iy,iz-dz) + F3_BASE(&psc.pf, HY,ix-dx,iy,iz-dz)))
#define HZ_CC(ix,iy,iz) (.25f*(F3_BASE(&psc.pf, HZ,ix,iy   ,iz) + F3_BASE(&psc.pf, HZ,ix-dx,iy   ,iz) + \
			       F3_BASE(&psc.pf, HZ,ix,iy-dy,iz) + F3_BASE(&psc.pf, HZ,ix-dx,iy-dy,iz)))

static void
calc_H(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = HX_CC(ix,iy,iz);
      F3_BASE(f, 1, ix,iy,iz) = HY_CC(ix,iy,iz);
      F3_BASE(f, 2, ix,iy,iz) = HZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

static void
calc_jdote(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
      F3_BASE(f, 1, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
      F3_BASE(f, 2, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

static void
calc_poyn(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = (EY_CC(ix,iy,iz) * HZ_CC(ix,iy,iz) - 
				 EZ_CC(ix,iy,iz) * HY_CC(ix,iy,iz));
      F3_BASE(f, 1, ix,iy,iz) = (EZ_CC(ix,iy,iz) * HX_CC(ix,iy,iz) -
				 EX_CC(ix,iy,iz) * HZ_CC(ix,iy,iz));
      F3_BASE(f, 2, ix,iy,iz) = (EX_CC(ix,iy,iz) * HY_CC(ix,iy,iz) -
				 EY_CC(ix,iy,iz) * HX_CC(ix,iy,iz));
    } foreach_3d_end;
  }
}

static void
calc_E2(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
      F3_BASE(f, 1, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
      F3_BASE(f, 2, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
}

static void
calc_H2(fields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(patch) {
    foreach_3d(patch, ix, iy, iz, 0, 0) {
      F3_BASE(f, 0, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
      F3_BASE(f, 1, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
      F3_BASE(f, 2, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
}

struct output_field {
  char *name;
  int nr_comp;
  char *fld_names[6];
  void (*calc)(fields_base_t *f);
};

static struct output_field output_fields[] = {
  { .name = "n"    , .nr_comp = 3, .fld_names = { "ne", "ni", "nn" },
    .calc = psc_calc_densities },
  { .name = "v"    , .nr_comp = 6, .fld_names = { "vex", "vey", "vez", "vix", "viy", "viz" },
    .calc = psc_calc_moments_v },
  { .name = "vv"    , .nr_comp = 6, .fld_names = { "vexvex", "veyvey", "vezvez",
						   "vixvix", "viyviy", "vizivz" },
    .calc = psc_calc_moments_vv },
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

    struct psc_patch *patch = &psc.patch[0];
    int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
    int ihg[3] = { patch->ldims[0] + psc.ibn[0],
		   patch->ldims[1] + psc.ibn[1],
		   patch->ldims[2] + psc.ibn[2] };
    fields_base_alloc(f, ilg, ihg, of->nr_comp);
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
    struct psc_patch *patch = &psc.patch[0];
    int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
    int ihg[3] = { patch->ldims[0] + psc.ibn[0],
		   patch->ldims[1] + psc.ibn[1],
		   patch->ldims[2] + psc.ibn[2] };
    fields_base_alloc(&tfd->flds[i], ilg, ihg, pfd->flds[i].nr_comp);
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
  &psc_output_format_ops_vtk_binary,
  &psc_output_format_ops_mrc,
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
  { "output_fields"      , VAR(output_fields)        , PARAM_STRING("n,j,e,h") },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)           },
  { "pfield_first"       , VAR(pfield_first)         , PARAM_INT(0)            },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)           },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)           },
  { "tfield_first"       , VAR(tfield_first)         , PARAM_INT(0)            },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)           },
  { "pfield_out_x_min"   , VAR(rn[0])                , PARAM_INT(0)            },  
  { "pfield_out_x_max"   , VAR(rx[0])                , PARAM_INT(1000000000)  },     // a big number to change it later to domain.ihi or command line number
  { "pfield_out_y_min"   , VAR(rn[1])                , PARAM_INT(0)           }, 
  { "pfield_out_y_max"   , VAR(rx[1])                , PARAM_INT(1000000000)  },
  { "pfield_out_z_min"   , VAR(rn[2])                , PARAM_INT(0)            }, 
  { "pfield_out_z_max"   , VAR(rx[2])                , PARAM_INT(1000000000)  },
  {},
};

#undef VAR

static struct psc_output_c psc_output_c;

// ----------------------------------------------------------------------
// output_c_create

static void output_c_create(void)
{ 
  struct psc_output_c *out = &psc_output_c;
  mrc_params_parse(out, psc_output_c_descr, "PSC output C", MPI_COMM_WORLD);
  mrc_params_print(out, psc_output_c_descr, "PSC output C", MPI_COMM_WORLD);

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
    struct psc_patch *patch = &psc.patch[0];
    for (int m = 0; m < f->nr_comp; m++) {
      fields_base_t *fld = &list->flds[list->nr_flds++];
      int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
      int ihg[3] = { patch->ldims[0] + psc.ibn[0],
		     patch->ldims[1] + psc.ibn[1],
		     patch->ldims[2] + psc.ibn[2] };
      fields_base_alloc_with_array(fld, ilg, ihg, 1,
				   &F3_BASE(f,m, -psc.ibn[0], -psc.ibn[1], -psc.ibn[2]));
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
// output_c_field

static void
output_c_field()
{
  struct psc_output_c *out = &psc_output_c;

  static bool first_time = true;
  struct psc_patch *patch = &psc.patch[0];
  if (first_time) {
    output_c_setup(out);
	  
	// set the output ranges
	for(int i=0;i<3;++i)
	{
	        if(out->rn[i]<0) out->rn[i]=0;
		if(out->rx[i]>psc.domain.gdims[i]) out->rx[i]=psc.domain.gdims[i];
		
		if(out->rx[i]>patch->off[i] + patch->ldims[i]) out->rx[i]=patch->off[i] + patch->ldims[i];
		if(out->rn[i]<patch->off[i]) out->rn[i]=patch->off[i];
		
		if(out->rn[i]>patch->off[i] + patch->ldims[i])
		{
		  out->rn[i]=patch->off[i] + patch->ldims[i];
		  out->rx[i]=out->rn[i];
		}
		if(out->rx[i]<patch->off[i]) 
		{
		  out->rx[i]=patch->off[i]; 
		  out->rn[i]=out->rx[i];
		}
		
	}
	
	// done setting output ranges
	  printf("rnx=%d\t rny=%d\t rnz=%d\n", out->rn[0], out->rn[1], out->rn[2]);
	  printf("rxx=%d\t rxy=%d\t rxz=%d\n", out->rx[0], out->rx[1], out->rx[2]);	  
	  
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
       out->format_ops->write_fields(out, &flds_list, "pfd");
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
      out->format_ops->write_fields(out, &flds_list, "tfd");
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

