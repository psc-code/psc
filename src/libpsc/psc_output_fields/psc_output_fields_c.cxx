
#include "psc_output_fields_c.h"
#include "psc_output_fields_item.h"
#include "psc_bnd.h"
#include "psc_fields_as_c.h"
#include "bnd.hxx"
#include "fields_item.hxx"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <string.h>

#define to_psc_output_fields_c(out) ((struct psc_output_fields_c *)((out)->obj.subctx))

void
write_fields(struct psc_output_fields_c *out, MfieldsList& list,
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
  mrc_domain_get_global_dims(ppsc->mrc_domain_, gdims);
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

  for (auto& item : list) {
    if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) {
      mrc_io_set_param_int3(io, "slab_off", slab_off);
      mrc_io_set_param_int3(io, "slab_dims", slab_dims);
    }

    auto mflds_base = item.mflds;
    auto& mflds = *item.mflds.sub();
    mflds.write_as_mrc_fld(io, item.name, item.comp_names);
  }
  mrc_io_close(io);
}

// ----------------------------------------------------------------------
// psc_output_fields_c_destroy

static void
psc_output_fields_c_destroy(struct psc_output_fields *out)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);

  for (int i = 0; i < out_c->item.size(); i++) {
    psc_output_fields_item_destroy(out_c->item[i]);
  }

  for (auto& item : out_c->tfd_) {
    psc_mfields_destroy(item.mflds.mflds());
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

  auto& pfd_ = out_c->pfd_, tfd_ = out_c->tfd_;
  
  // setup pfd according to output_fields as given
  // (potentially) on the command line
  // parse comma separated list of fields
  char *s_orig = strdup(out_c->output_fields), *p, *s = s_orig;
  while ((p = strsep(&s, ", "))) {
    struct psc_output_fields_item *item =
      psc_output_fields_item_create(psc_output_fields_comm(out));
    psc_output_fields_item_set_type(item, p);
    psc_output_fields_item_setup(item);
    out_c->item.push_back(item);

    // pfd_
    auto mres = PscFieldsItemBase{item}->mres();
    std::vector<std::string> comp_names;
    for (int m = 0; m < mres->n_comps(); m++) {
      comp_names.push_back(psc_mfields_comp_name(mres.mflds(), m));
    }
    pfd_.emplace_back(mres, p, comp_names);

    // tfd_ -- FIXME?! always MfieldsC
    auto mflds = PscMfields<MfieldsC>::create(psc_comm(psc), psc->grid(), mres->n_comps(), psc->ibn);
    tfd_.emplace_back(PscMfieldsBase{mflds.mflds()}, p, comp_names);
  }
  free(s_orig);

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

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			struct psc_mfields *flds, struct psc_mparticles *particles)
{
  struct psc_output_fields_c *out_c = to_psc_output_fields_c(out);
  PscMparticlesBase mprts(particles);
  struct psc *psc = out->psc;

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  bool doaccum_tfield = out_c->dowrite_tfield && 
    (((psc->timestep >= (out_c->tfield_next - out_c->tfield_length + 1)) &&
      psc->timestep % out_c->tfield_every == 0) ||
     psc->timestep == 0);

  auto& pfd_ = out_c->pfd_, tfd_ = out_c->tfd_;
  
  if ((out_c->dowrite_pfield && psc->timestep >= out_c->pfield_next) ||
      (out_c->dowrite_tfield && doaccum_tfield)) {
    for (int i = 0; i < out_c->item.size(); i++) {
      PscFieldsItemBase item(out_c->item[i]);
      item(flds, mprts);
    }
  }

  if (out_c->dowrite_pfield && psc->timestep >= out_c->pfield_next) {
    mpi_printf(psc_output_fields_comm(out), "***** Writing PFD output\n");
    out_c->pfield_next += out_c->pfield_step;
    write_fields(out_c, pfd_, IO_TYPE_PFD, out_c->pfd_s);
  }

  if (out_c->dowrite_tfield) {
    if (doaccum_tfield) {
      // tfd += pfd
      for (int m = 0; m < tfd_.size(); m++) {
	tfd_[m].mflds->axpy(1., *pfd_[m].mflds.sub());
      }
      out_c->naccum++;
    }
    if (psc->timestep >= out_c->tfield_next) {
      mpi_printf(psc_output_fields_comm(out), "***** Writing TFD output\n");
      out_c->tfield_next += out_c->tfield_step;
      
      // convert accumulated values to correct temporal mean
      for (auto item: tfd_) {
	item.mflds->scale(1. / out_c->naccum);
      }
      
      write_fields(out_c, tfd_, IO_TYPE_TFD, out_c->tfd_s);
      
      for (auto item: tfd_) {
	item.mflds->zero();
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
  { "output_fields"      , VAR(output_fields)        , PARAM_STRING("j,e,h")   },
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

struct psc_output_fields_ops_c : psc_output_fields_ops {
  psc_output_fields_ops_c() {
    name                  = "c";
    size                  = sizeof(struct psc_output_fields_c);
    param_descr           = psc_output_fields_c_descr;
    setup                 = psc_output_fields_c_setup;
    destroy               = psc_output_fields_c_destroy;
    write                 = psc_output_fields_c_write;
    read                  = psc_output_fields_c_read;
    run                   = psc_output_fields_c_run;
  }
} psc_output_fields_c_ops;
