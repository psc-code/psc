
#include "psc_output_fields_c.h"
#include "psc_output_fields_item.h"
#include "psc_bnd.h"
#include "psc_fields_as_c.h"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <string.h>

#define to_psc_output_fields_c(out) ((struct psc_output_fields_c *)((out)->obj.subctx))

void
write_fields(struct psc_output_fields_c *out, struct psc_fields_list *list,
	     int io_type, const char *pfx)
{
  struct mrc_io *io = out->ios[io_type];
  if (!io) {
    io = mrc_io_create(MPI_COMM_WORLD);
    mrc_io_set_param_string(io, "basename", pfx);
    mrc_io_set_param_string(io, "outdir",out->data_dir);
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
    mrc_io_view(io);
    out->ios[io_type] = io;
  }

  int gdims[3];
  mrc_domain_get_global_dims(ppsc->mrc_domain, gdims);
  int slab_off[3], slab_dims[3];
  for (int d = 0; d < 3; d++) {
    if (out->rx[d] > gdims[d])
      out->rx[d] = gdims[d];
    
    slab_off[d] = out->rn[d];
    slab_dims[d] = out->rx[d] - out->rn[d];
  }

  mrc_io_open(io, "w", ppsc->timestep, ppsc->timestep * ppsc->dt);

  // save some basic info about the run in the output file
  struct mrc_obj *obj = mrc_obj_create(mrc_io_comm(io));
  mrc_obj_set_name(obj, "psc");
  mrc_obj_dict_add_int(obj, "timestep", ppsc->timestep);
  mrc_obj_dict_add_float(obj, "time", ppsc->timestep * ppsc->dt);
  mrc_obj_dict_add_float(obj, "cc", ppsc->prm.cc);
  mrc_obj_dict_add_float(obj, "dt", ppsc->dt);
  mrc_obj_write(obj, io);
  mrc_obj_destroy(obj);

  for (int m = 0; m < list->nr_flds; m++) {
    struct psc_mfields *flds = list->flds[m];

    if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) {
      mrc_io_set_param_int3(io, "slab_off", slab_off);
      mrc_io_set_param_int3(io, "slab_dims", slab_dims);
    }

    psc_mfields_write_as_mrc_fld(flds, io);
  }
  mrc_io_close(io);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_create

static void
psc_output_fields_c_create(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

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

  // FIXME, we used to have mrc_io persistent across new objects from
  // this class (using a static var for ios), which was bad, but actually
  // helped getting a proper temporal XDMF file...
  // now it's part of the object, but will break the temporal XDMF file on
  // rebalancing (?)
  
  for (int i = 0; i < NR_IO_TYPES; i++) {
    mrc_io_destroy(out_c->ios[i]);
    out_c->ios[i] = NULL;
  }
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
    psc_output_fields_item_setup(item);
    out_c->item[pfd->nr_flds] = item;
    struct psc_mfields *flds = psc_output_fields_item_create_mfields(item);
    psc_mfields_set_name(flds, p);
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
    // FIXME, shouldn't we use item_create_mfields(), too?
    struct psc_mfields *flds = psc_mfields_create(mrc_domain_comm(psc->mrc_domain));
    psc_mfields_set_type(flds, "c");
    psc_mfields_set_name(flds, psc_mfields_name(pfd->flds[i]));
    psc_mfields_set_param_obj(flds, "domain", psc->mrc_domain);
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
// psc_output_fields_c_write

static void
psc_output_fields_c_write(struct psc_output_fields *out, struct mrc_io *io)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  mrc_io_write_int(io, out, "pfield_next", out_c->pfield_next);
  mrc_io_write_int(io, out, "tfield_next", out_c->tfield_next);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_read

static void
psc_output_fields_c_read(struct psc_output_fields *out, struct mrc_io *io)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  psc_output_fields_read_super(out, io);
  mrc_io_read_int(io, out, "pfield_next", &out_c->pfield_next);
  mrc_io_read_int(io, out, "tfield_next", &out_c->tfield_next);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_run

bool psc_output_fields_check_bnd;
void psc_bnd_check_domain(struct psc_bnd *bnd);

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			struct psc_mfields *flds, struct psc_mparticles *particles)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  struct psc *psc = out->psc;

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  if (psc_output_fields_check_bnd) {
    psc_output_fields_check_bnd = false;
    psc_bnd_check_domain(out_c->bnd);
  }

  bool doaccum_tfield = out_c->dowrite_tfield && 
    (((psc->timestep >= (out_c->tfield_next - out_c->tfield_length + 1)) &&
      psc->timestep % out_c->tfield_every == 0) || 
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
       write_fields(out_c, &out_c->pfd, IO_TYPE_PFD, out_c->pfd_s);
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

      write_fields(out_c, &out_c->tfd, IO_TYPE_TFD, out_c->tfd_s);

      for (int m = 0; m < out_c->tfd.nr_flds; m++) {
	for (int mm = 0; mm < out_c->tfd.flds[m]->nr_fields; mm++) {
	  psc_mfields_zero_comp(out_c->tfd.flds[m], mm);
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
  { "pfd"                , VAR(pfd_s)                , PARAM_STRING("pfd")     },
  { "tfd"                , VAR(tfd_s)                , PARAM_STRING("tfd")     },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)           },
  { "pfield_first"       , VAR(pfield_first)         , PARAM_INT(0)            },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)           },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)           },
  { "tfield_first"       , VAR(tfield_first)         , PARAM_INT(0)            },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)           },
  { "tfield_length"      , VAR(tfield_length)        , PARAM_INT(10)           },
  { "tfield_every"       , VAR(tfield_every)         , PARAM_INT(1)            },
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
  .destroy               = psc_output_fields_c_destroy,
  .write                 = psc_output_fields_c_write,
  .read                  = psc_output_fields_c_read,
  .run                   = psc_output_fields_c_run,
};
