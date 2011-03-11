
#include "psc_output_fields_c.h"

#include <mrc_profile.h>

#define to_psc_output_fields_c(out) ((struct psc_output_fields_c *)((out)->obj.subctx))

// ======================================================================

#define define_dxdydz(dx, dy, dz)					\
  int dx __unused = (psc.domain.gdims[0] == 1) ? 0 : 1;			\
  int dy __unused = (psc.domain.gdims[1] == 1) ? 0 : 1;			\
  int dz __unused = (psc.domain.gdims[2] == 1) ? 0 : 1

#define JX_CC(ix,iy,iz) (.5f * (F3_BASE(pf, JXI,ix,iy,iz) + F3_BASE(pf, JXI,ix-dx,iy,iz)))
#define JY_CC(ix,iy,iz) (.5f * (F3_BASE(pf, JYI,ix,iy,iz) + F3_BASE(pf, JYI,ix,iy-dy,iz)))
#define JZ_CC(ix,iy,iz) (.5f * (F3_BASE(pf, JZI,ix,iy,iz) + F3_BASE(pf, JZI,ix,iy,iz-dz)))

static void
calc_j(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = JX_CC(ix,iy,iz);
      F3_BASE(ff, 1, ix,iy,iz) = JY_CC(ix,iy,iz);
      F3_BASE(ff, 2, ix,iy,iz) = JZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

#define EX_CC(ix,iy,iz) (.5f * (F3_BASE(pf, EX,ix,iy,iz) + F3_BASE(pf, EX,ix-dx,iy,iz)))
#define EY_CC(ix,iy,iz) (.5f * (F3_BASE(pf, EY,ix,iy,iz) + F3_BASE(pf, EY,ix,iy-dy,iz)))
#define EZ_CC(ix,iy,iz) (.5f * (F3_BASE(pf, EZ,ix,iy,iz) + F3_BASE(pf, EZ,ix,iy,iz-dz)))

static void
calc_E(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = EX_CC(ix,iy,iz);
      F3_BASE(ff, 1, ix,iy,iz) = EY_CC(ix,iy,iz);
      F3_BASE(ff, 2, ix,iy,iz) = EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

#define HX_CC(ix,iy,iz) (.25f*(F3_BASE(pf, HX,ix,iy,iz   ) + F3_BASE(pf, HX,ix,iy-dy,iz   ) + \
			       F3_BASE(pf, HX,ix,iy,iz-dz) + F3_BASE(pf, HX,ix,iy-dy,iz-dz)))
#define HY_CC(ix,iy,iz) (.25f*(F3_BASE(pf, HY,ix,iy,iz   ) + F3_BASE(pf, HY,ix-dx,iy,iz   ) + \
			       F3_BASE(pf, HY,ix,iy,iz-dz) + F3_BASE(pf, HY,ix-dx,iy,iz-dz)))
#define HZ_CC(ix,iy,iz) (.25f*(F3_BASE(pf, HZ,ix,iy   ,iz) + F3_BASE(pf, HZ,ix-dx,iy   ,iz) + \
			       F3_BASE(pf, HZ,ix,iy-dy,iz) + F3_BASE(pf, HZ,ix-dx,iy-dy,iz)))

static void
calc_H(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = HX_CC(ix,iy,iz);
      F3_BASE(ff, 1, ix,iy,iz) = HY_CC(ix,iy,iz);
      F3_BASE(ff, 2, ix,iy,iz) = HZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

static void
calc_jdote(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = JX_CC(ix,iy,iz) * EX_CC(ix,iy,iz);
      F3_BASE(ff, 1, ix,iy,iz) = JY_CC(ix,iy,iz) * EY_CC(ix,iy,iz);
      F3_BASE(ff, 2, ix,iy,iz) = JZ_CC(ix,iy,iz) * EZ_CC(ix,iy,iz);
    } foreach_3d_end;
  }
}

static void
calc_poyn(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = (EY_CC(ix,iy,iz) * HZ_CC(ix,iy,iz) - 
				  EZ_CC(ix,iy,iz) * HY_CC(ix,iy,iz));
      F3_BASE(ff, 1, ix,iy,iz) = (EZ_CC(ix,iy,iz) * HX_CC(ix,iy,iz) -
				  EX_CC(ix,iy,iz) * HZ_CC(ix,iy,iz));
      F3_BASE(ff, 2, ix,iy,iz) = (EX_CC(ix,iy,iz) * HY_CC(ix,iy,iz) -
				  EY_CC(ix,iy,iz) * HX_CC(ix,iy,iz));
    } foreach_3d_end;
  }
}

static void
calc_E2(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = sqr(EX_CC(ix,iy,iz));
      F3_BASE(ff, 1, ix,iy,iz) = sqr(EY_CC(ix,iy,iz));
      F3_BASE(ff, 2, ix,iy,iz) = sqr(EZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
}

static void
calc_H2(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f)
{
  define_dxdydz(dx, dy, dz);
  foreach_patch(p) {
    fields_base_t *ff = &f->f[p];
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 0, 0) {
      F3_BASE(ff, 0, ix,iy,iz) = sqr(HX_CC(ix,iy,iz));
      F3_BASE(ff, 1, ix,iy,iz) = sqr(HY_CC(ix,iy,iz));
      F3_BASE(ff, 2, ix,iy,iz) = sqr(HZ_CC(ix,iy,iz));
    } foreach_3d_end;
  }
}

struct output_field {
  char *name;
  int nr_comp;
  char *fld_names[6];
  void (*calc)(mfields_base_t *flds, mparticles_base_t *particles, mfields_base_t *f);
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

// ======================================================================
// find_output_format_ops

// FIXME, this should be converted to mrc_obj, too

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

// ----------------------------------------------------------------------
// psc_output_fields_c_create

static void
psc_output_fields_c_create(struct mrc_obj *obj)
{
  struct psc_output_fields *out = to_psc_output_fields(obj);
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  out_c->pfield_next = out_c->pfield_first;
  out_c->tfield_next = out_c->tfield_first;

  out_c->format_ops = find_output_format_ops(out_c->output_format);
  if (out_c->format_ops->create) {
    out_c->format_ops->create();
  }
}

// ----------------------------------------------------------------------
// output_c_setup

static void
output_c_setup(struct psc_output_fields_c *out)
{
  struct psc_fields_list *pfd = &out->pfd;

  // setup pfd according to output_fields as given
  // (potentially) on the command line
  pfd->nr_flds = 0;
  // parse comma separated list of fields
  char *p, *s = strdup(out->output_fields);
  while ((p = strsep(&s, ", "))) {
    struct output_field *of = find_output_field(p);
    mfields_base_t *flds = &pfd->flds[pfd->nr_flds];
    out->out_flds[pfd->nr_flds] = of;
    pfd->nr_flds++;

    mfields_base_alloc(flds, of->nr_comp);
    foreach_patch(pp) {
      for (int m = 0; m < of->nr_comp; m++) {
	flds->f[pp].name[m] = strdup(of->fld_names[m]);
      }
    }
  }
  free(s);

  // create tfd to look just like pfd
  // FIXME, only if necessary
  struct psc_fields_list *tfd = &out->tfd;
  tfd->nr_flds = pfd->nr_flds;
  for (int i = 0; i < pfd->nr_flds; i++) {
    mfields_base_alloc(&tfd->flds[i], pfd->flds[i].f[0].nr_comp);
    foreach_patch(pp) {
      for (int m = 0; m < pfd->flds[i].f[pp].nr_comp; m++) {
	tfd->flds[i].f[pp].name[m] = strdup(pfd->flds[i].f[pp].name[m]);
      }
    }
  }
  out->naccum = 0;
}

// ----------------------------------------------------------------------
// make_fields_list

static void
make_fields_list(struct psc_fields_list *list, struct psc_fields_list *list_in)
{
  // the only thing this still does is to flatten
  // the list so that it only contains 1-component entries

  list->nr_flds = 0;
  for (int i = 0; i < list_in->nr_flds; i++) {
    mfields_base_t *flds_in = &list_in->flds[i];
    for (int m = 0; m < flds_in->f[0].nr_comp; m++) {
      mfields_base_t *flds = &list->flds[list->nr_flds++];
      flds->f = calloc(psc.nr_patches, sizeof(*flds->f));
      foreach_patch(p) {
	int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
	int ihg[3] = { psc.patch[p].ldims[0] + psc.ibn[0],
		       psc.patch[p].ldims[1] + psc.ibn[1],
		       psc.patch[p].ldims[2] + psc.ibn[2] };
	fields_base_alloc_with_array(&flds->f[p], ilg, ihg, 1,
				     &F3_BASE(&flds_in->f[p],m, -psc.ibn[0], -psc.ibn[1], -psc.ibn[2]));
	flds->f[p].name[0] = strdup(flds_in->f[p].name[m]);
      }
    }
  }
}

// ----------------------------------------------------------------------
// free_fields_list

static void
free_fields_list(struct psc_fields_list *list)
{
  for (int m = 0; m < list->nr_flds; m++) {
    foreach_patch(p) {
      fields_base_free(&list->flds[m].f[p]);
    }
    free(list->flds[m].f);
  }
}

// ----------------------------------------------------------------------
// psc_output_fields_c_run

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			mfields_base_t *flds, mparticles_base_t *particles)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  static bool first_time = true;
  struct psc_patch *patch = &psc.patch[0];

  if (first_time) {
    output_c_setup(out_c);
	  
    // set the output ranges
    for(int i=0;i<3;++i) {
      if(out_c->rn[i]<0) out_c->rn[i]=0;
      if(out_c->rx[i]>psc.domain.gdims[i]) out_c->rx[i]=psc.domain.gdims[i];
      
      if(out_c->rx[i]>patch->off[i] + patch->ldims[i]) out_c->rx[i]=patch->off[i] + patch->ldims[i];
      if(out_c->rn[i]<patch->off[i]) out_c->rn[i]=patch->off[i];
      
      if(out_c->rn[i]>patch->off[i] + patch->ldims[i]) {
	out_c->rn[i]=patch->off[i] + patch->ldims[i];
	out_c->rx[i]=out_c->rn[i];
      }
      if(out_c->rx[i]<patch->off[i]) {
	out_c->rx[i]=patch->off[i]; 
	out_c->rn[i]=out_c->rx[i];
      }
    }
	
    // done setting output ranges
    printf("rnx=%d\t rny=%d\t rnz=%d\n", out_c->rn[0], out_c->rn[1], out_c->rn[2]);
    printf("rxx=%d\t rxy=%d\t rxz=%d\n", out_c->rx[0], out_c->rx[1], out_c->rx[2]);	  
	  
    first_time = false;
  }

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  if ((out_c->dowrite_pfield && psc.timestep >= out_c->pfield_next) ||
      out_c->dowrite_tfield) {
    struct psc_fields_list *pfd = &out_c->pfd;
    for (int i = 0; i < pfd->nr_flds; i++) {
      out_c->out_flds[i]->calc(flds, particles, &pfd->flds[i]);
    }
  }
  
  if (out_c->dowrite_pfield) {
    if (psc.timestep >= out_c->pfield_next) {
       out_c->pfield_next += out_c->pfield_step;
       struct psc_fields_list flds_list;
       make_fields_list(&flds_list, &out_c->pfd);
       out_c->format_ops->write_fields(out_c, &flds_list, "pfd");
       free_fields_list(&flds_list);
    }
  }

  if (out_c->dowrite_tfield) {
    foreach_patch(p) {
      for (int m = 0; m < out_c->tfd.nr_flds; m++) {
	// tfd += pfd
	fields_base_axpy_all(&out_c->tfd.flds[m].f[p], 1., &out_c->pfd.flds[m].f[p]);
      }
    }
    out_c->naccum++;
    if (psc.timestep >= out_c->tfield_next) {
      out_c->tfield_next += out_c->tfield_step;

      // convert accumulated values to correct temporal mean
      foreach_patch(p) {
	for (int m = 0; m < out_c->tfd.nr_flds; m++) {
	  fields_base_scale_all(&out_c->tfd.flds[m].f[p], 1. / out_c->naccum);
	}
      }

      struct psc_fields_list flds_list;
      make_fields_list(&flds_list, &out_c->tfd);
      out_c->format_ops->write_fields(out_c, &flds_list, "tfd");
      free_fields_list(&flds_list);
      foreach_patch(p) {
	for (int m = 0; m < out_c->tfd.nr_flds; m++) {
	  fields_base_zero_all(&out_c->tfd.flds[m].f[p]);
	}
      }
      out_c->naccum = 0;
    }
  }
  
  prof_stop(pr);
}

// ======================================================================
// psc_output_fields: subclass "c"

#define VAR(x) (void *)offsetof(struct psc_output_fields_c, x)

static struct param psc_output_fields_c_descr[] = {
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

struct psc_output_fields_ops psc_output_fields_c_ops = {
  .name                  = "c",
  .size                  = sizeof(struct psc_output_fields_c),
  .param_descr           = psc_output_fields_c_descr,
  .create                = psc_output_fields_c_create,
  .run                   = psc_output_fields_c_run,
};
