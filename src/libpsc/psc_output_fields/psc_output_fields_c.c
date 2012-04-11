
#include "psc_output_fields_c.h"
#include "psc_output_fields_item.h"
#include "psc_output_format.h"
#include "psc_bnd.h"
#include "psc_fields_as_c.h"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <string.h>

#define to_psc_output_fields_c(out) ((struct psc_output_fields_c *)((out)->obj.subctx))

// ----------------------------------------------------------------------
// psc_output_fields_c_create

static void
psc_output_fields_c_create(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  out_c->format = psc_output_format_create(psc_output_fields_comm(out));

  out_c->bnd = psc_bnd_create(psc_output_fields_comm(out));
  psc_bnd_set_name(out_c->bnd, "psc_output_fields_bnd");
  psc_bnd_set_type(out_c->bnd, "c");
  psc_bnd_set_psc(out_c->bnd, ppsc);
  psc_output_fields_add_child(out, (struct mrc_obj *) out_c->bnd);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_destroy

static void
psc_output_fields_c_destroy(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  struct psc_fields_list *pfd = &out_c->pfd;
  for (int i = 0; i < pfd->nr_flds; i++) {
    psc_mfields_list_del(&psc_mfields_base_list, &pfd->flds[i]);
    psc_mfields_destroy(pfd->flds[i]);
    psc_output_fields_item_destroy(out_c->item[i]);
  }
  struct psc_fields_list *tfd = &out_c->tfd;
  for (int i = 0; i < tfd->nr_flds; i++) {
    psc_mfields_list_del(&psc_mfields_base_list, &tfd->flds[i]);
    psc_mfields_destroy(tfd->flds[i]);
  }

  psc_output_format_destroy(out_c->format);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_set_from_options

static void
psc_output_fields_c_set_from_options(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  psc_output_format_set_from_options(out_c->format);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_setup

static void
psc_output_fields_c_setup(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  struct psc *psc = out->psc;

  out_c->pfield_next = out_c->pfield_first;
  out_c->tfield_next = out_c->tfield_first;

  struct psc_fields_list *pfd = &out_c->pfd;

  // setup pfd according to output_fields as given
  // (potentially) on the command line
  pfd->nr_flds = 0;
  // parse comma separated list of fields
  char *s_orig = strdup(out_c->output_fields), *p, *s = s_orig;
  while ((p = strsep(&s, ", "))) {
    struct psc_output_fields_item *item =
      psc_output_fields_item_create(psc_output_fields_comm(out));
    psc_output_fields_item_set_type(item, p);
    psc_output_fields_item_set_psc_bnd(item, out_c->bnd);
    out_c->item[pfd->nr_flds] = item;
    mfields_c_t *flds = psc_output_fields_item_create_mfields(item);
    pfd->flds[pfd->nr_flds] = flds;
    // FIXME, should be del'd eventually
    psc_mfields_list_add(&psc_mfields_base_list, &pfd->flds[pfd->nr_flds]);
    pfd->nr_flds++;
  }
  free(s_orig);

  // create tfd to look just like pfd
  // FIXME, only if necessary
  struct psc_fields_list *tfd = &out_c->tfd;
  tfd->nr_flds = pfd->nr_flds;
  for (int i = 0; i < pfd->nr_flds; i++) {
    assert(psc->nr_patches > 0);
    mfields_c_t *flds = psc_mfields_create(mrc_domain_comm(psc->mrc_domain));
    psc_mfields_set_type(flds, "c");
    psc_mfields_set_domain(flds, psc->mrc_domain);
    psc_mfields_set_param_int(flds, "nr_fields", pfd->flds[i]->nr_fields);
    psc_mfields_set_param_int3(flds, "ibn", psc->ibn);
    psc_mfields_setup(flds);
    tfd->flds[i] = flds;
    // FIXME, should be del'd eventually
    psc_mfields_list_add(&psc_mfields_base_list, &tfd->flds[i]);
    for (int m = 0; m < pfd->flds[i]->nr_fields; m++) {
      psc_mfields_set_comp_name(flds, m, psc_mfields_comp_name(pfd->flds[i], m));
    }
  }
  out_c->naccum = 0;

  psc_output_fields_setup_children(out);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_view

static void
psc_output_fields_c_view(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  psc_output_format_view(out_c->format);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_write

static void
psc_output_fields_c_write(struct psc_output_fields *out, struct mrc_io *io)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  const char *path = psc_output_fields_name(out);
  mrc_io_write_attr_int(io, path, "pfield_next", out_c->pfield_next);
  mrc_io_write_attr_int(io, path, "tfield_next", out_c->tfield_next);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_read

static void
psc_output_fields_c_read(struct psc_output_fields *out, struct mrc_io *io)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  psc_output_fields_read_super(out, io);

  const char *path = psc_output_fields_name(out);
  mrc_io_read_attr_int(io, path, "pfield_next", &out_c->pfield_next);
  mrc_io_read_attr_int(io, path, "tfield_next", &out_c->tfield_next);
}

// ----------------------------------------------------------------------
// make_fields_list

static void
make_fields_list(struct psc *psc, struct psc_fields_list *list,
		 struct psc_fields_list *list_in)
{
  // the only thing this still does is to flatten
  // the list so that it only contains 1-component entries
  // FIXME, slow and unnec

  list->nr_flds = 0;
  for (int i = 0; i < list_in->nr_flds; i++) {
    mfields_c_t *flds_in = list_in->flds[i];
    for (int m = 0; m < flds_in->nr_fields; m++) {
      mfields_c_t *flds = psc_mfields_create(psc_comm(psc));
      psc_mfields_set_type(flds, "c");
      psc_mfields_set_domain(flds, psc->mrc_domain);
      psc_mfields_set_param_int3(flds, "ibn", psc->ibn);
      psc_mfields_setup(flds);
      psc_mfields_copy_comp(flds, 0, flds_in, m);
      list->flds[list->nr_flds++] = flds;
      psc_mfields_set_comp_name(flds, 0, psc_mfields_comp_name(flds_in, m));
    }
  }
}

// ----------------------------------------------------------------------
// free_fields_list

static void
free_fields_list(struct psc *psc, struct psc_fields_list *list)
{
  for (int m = 0; m < list->nr_flds; m++) {
    psc_mfields_destroy(list->flds[m]);
  }
}

// ----------------------------------------------------------------------
// psc_output_fields_c_run

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			mfields_base_t *flds, mparticles_base_t *particles)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  struct psc *psc = out->psc;

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  bool doaccum_tfield = out_c->dowrite_tfield && 
        ((psc->timestep >= (out_c->tfield_next - out_c->tfield_length + 1)) || 
          psc->timestep == 0);

  if ((out_c->dowrite_pfield && psc->timestep >= out_c->pfield_next) ||
      (out_c->dowrite_tfield && doaccum_tfield)) {
    struct psc_fields_list *pfd = &out_c->pfd;
    for (int i = 0; i < pfd->nr_flds; i++) {
      psc_output_fields_item_run(out_c->item[i], flds, particles, pfd->flds[i]);
    }
  }
  
  if (out_c->dowrite_pfield) {
    if (psc->timestep >= out_c->pfield_next) {
       out_c->pfield_next += out_c->pfield_step;
       struct psc_fields_list flds_list;
       make_fields_list(psc, &flds_list, &out_c->pfd);
       psc_output_format_write_fields(out_c->format, out_c, &flds_list, "pfd");
       free_fields_list(psc, &flds_list);
    }
  }

  if (out_c->dowrite_tfield) {
   if (doaccum_tfield) {
    // tfd += pfd
    for (int m = 0; m < out_c->tfd.nr_flds; m++) {
      psc_mfields_axpy(out_c->tfd.flds[m], 1., out_c->pfd.flds[m]);
    }
    out_c->naccum++;
   }
    if (psc->timestep >= out_c->tfield_next) {
      out_c->tfield_next += out_c->tfield_step;

      // convert accumulated values to correct temporal mean
      for (int m = 0; m < out_c->tfd.nr_flds; m++) {
	psc_mfields_scale(out_c->tfd.flds[m], 1. / out_c->naccum);
      }

      struct psc_fields_list flds_list;
      make_fields_list(psc, &flds_list, &out_c->tfd);
      psc_output_format_write_fields(out_c->format, out_c, &flds_list, "tfd");
      free_fields_list(psc, &flds_list);
      for (int m = 0; m < out_c->tfd.nr_flds; m++) {
	for (int mm = 0; mm < out_c->tfd.flds[m]->nr_fields; mm++) {
	  psc_mfields_zero(out_c->tfd.flds[m], mm);
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

// FIXME pfield_out_[xyz]_{min,max} aren't for pfield only, better init to 0,
// use INT3

static struct param psc_output_fields_c_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "output_fields"      , VAR(output_fields)        , PARAM_STRING("n,j,e,h") },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)           },
  { "pfield_first"       , VAR(pfield_first)         , PARAM_INT(0)            },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)           },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)           },
  { "tfield_first"       , VAR(tfield_first)         , PARAM_INT(0)            },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)           },
  { "tfield_length"      , VAR(tfield_length)        , PARAM_INT(10)           },
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
  .setup                 = psc_output_fields_c_setup,
  .set_from_options      = psc_output_fields_c_set_from_options,
  .destroy               = psc_output_fields_c_destroy,
  .write                 = psc_output_fields_c_write,
  .read                  = psc_output_fields_c_read,
  .view                  = psc_output_fields_c_view,
  .run                   = psc_output_fields_c_run,
};
