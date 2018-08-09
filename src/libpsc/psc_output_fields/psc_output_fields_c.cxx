
#include "psc_output_fields_c.h"
#include "psc_output_fields_item.h"
#include "psc_fields_as_c.h"
#include "bnd.hxx"
#include "fields_item.hxx"
#include "output_fields.hxx"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <string.h>

using PscOutputFields_t = PscOutputFields<psc_output_fields_c>;

static mrc_io* create_mrc_io(const char* pfx, const char* data_dir)
{
  mrc_io* io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_param_string(io, "basename", pfx);
  mrc_io_set_param_string(io, "outdir", data_dir);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);
  return io;
}

static void open_mrc_io(PscOutputFields_t outf, mrc_io *io)
{
  int gdims[3];
  mrc_domain_get_global_dims(ppsc->mrc_domain_, gdims);
  int slab_off[3], slab_dims[3];
  for (int d = 0; d < 3; d++) {
    if (outf->rx[d] > gdims[d])
      outf->rx[d] = gdims[d];
    
    slab_off[d] = outf->rn[d];
    slab_dims[d] = outf->rx[d] - outf->rn[d];
  }

  mrc_io_open(io, "w", ppsc->timestep, ppsc->timestep * ppsc->grid().dt);

  // save some basic info about the run in the output file
  struct mrc_obj *obj = mrc_obj_create(mrc_io_comm(io));
  mrc_obj_set_name(obj, "psc");
  mrc_obj_dict_add_int(obj, "timestep", ppsc->timestep);
  mrc_obj_dict_add_float(obj, "time", ppsc->timestep * ppsc->grid().dt);
  mrc_obj_dict_add_float(obj, "cc", ppsc->grid().norm.cc);
  mrc_obj_dict_add_float(obj, "dt", ppsc->grid().dt);
  mrc_obj_write(obj, io);
  mrc_obj_destroy(obj);

  if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) {
    mrc_io_set_param_int3(io, "slab_off", slab_off);
    mrc_io_set_param_int3(io, "slab_dims", slab_dims);
  }
}

// ----------------------------------------------------------------------
// psc_output_fields_c_destroy

static void
psc_output_fields_c_destroy(struct psc_output_fields *out)
{
  PscOutputFields_t outf{out};

  for (auto& item : outf->items) {
    psc_output_fields_item_destroy(item.item.item());
    delete &item.tfd;
  }

  for (int i = 0; i < NR_IO_TYPES; i++) {
    mrc_io_destroy(outf->ios[i]);
    outf->ios[i] = NULL;
  }
}

// ----------------------------------------------------------------------
// psc_output_fields_c_setup

static void
psc_output_fields_c_setup(struct psc_output_fields *out)
{
  PscOutputFields_t outf{out};
  struct psc *psc = ppsc;

  outf->pfield_next = outf->pfield_first;
  outf->tfield_next = outf->tfield_first;

  // setup pfd according to output_fields as given
  // (potentially) on the command line
  // parse comma separated list of fields
  char *s_orig = strdup(outf->output_fields), *p, *s = s_orig;
  while ((p = strsep(&s, ", "))) {
    struct psc_output_fields_item *item =
      psc_output_fields_item_create(psc_output_fields_comm(out));
    psc_output_fields_item_set_type(item, p);
    psc_output_fields_item_setup(item);

    // pfd
    std::vector<std::string> comp_names = PscFieldsItemBase{item}->comp_names();
    MfieldsBase& mflds_pfd = PscFieldsItemBase{item}->mres();

    // tfd -- FIXME?! always MfieldsC
    MfieldsBase& mflds_tfd = *new MfieldsC{psc->grid(), mflds_pfd.n_comps(), psc->ibn};
    outf->items.emplace_back(PscFieldsItemBase{item}, p, comp_names, mflds_pfd, mflds_tfd);
  }
  free(s_orig);

  outf->naccum = 0;

  psc_output_fields_setup_children(out);

  if (outf->dowrite_pfield) {
    outf->ios[IO_TYPE_PFD] = create_mrc_io(outf->pfd_s, outf->data_dir);
  }
  if (outf->dowrite_tfield) {
    outf->ios[IO_TYPE_TFD] = create_mrc_io(outf->tfd_s, outf->data_dir);
  }
}

// ----------------------------------------------------------------------
// psc_output_fields_c_run

static void
psc_output_fields_c_run(struct psc_output_fields *out,
			MfieldsBase& mflds, MparticlesBase& mprts)
{
  PscOutputFields_t outf{out};
  struct psc *psc = ppsc;

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  bool doaccum_tfield = outf->dowrite_tfield && 
    (((psc->timestep >= (outf->tfield_next - outf->tfield_length + 1)) &&
      psc->timestep % outf->tfield_every == 0) ||
     psc->timestep == 0);

  if ((outf->dowrite_pfield && psc->timestep >= outf->pfield_next) ||
      (outf->dowrite_tfield && doaccum_tfield)) {
    for (auto item : outf->items) {
      item.item(mflds, mprts);
    }
  }

  if (outf->dowrite_pfield && psc->timestep >= outf->pfield_next) {
    mpi_printf(psc_output_fields_comm(out), "***** Writing PFD output\n");
    outf->pfield_next += outf->pfield_step;

    auto io = outf->ios[IO_TYPE_PFD];
    open_mrc_io(outf, io);
    for (auto& item : outf->items) {
      item.pfd.write_as_mrc_fld(io, item.name, item.comp_names);
    }
    mrc_io_close(io);
  }

  if (outf->dowrite_tfield) {
    if (doaccum_tfield) {
      // tfd += pfd
      for (auto& item : outf->items) {
	item.tfd.axpy(1., item.pfd);
      }
      outf->naccum++;
    }
    if (psc->timestep >= outf->tfield_next) {
      mpi_printf(psc_output_fields_comm(out), "***** Writing TFD output\n");
      outf->tfield_next += outf->tfield_step;
      
      auto io = outf->ios[IO_TYPE_TFD];
      open_mrc_io(outf, io);

      // convert accumulated values to correct temporal mean
      for (auto& item : outf->items) {
	item.tfd.scale(1. / outf->naccum);
	item.tfd.write_as_mrc_fld(io, item.name, item.comp_names);
	item.tfd.zero();
      }
      outf->naccum = 0;
      mrc_io_close(io);
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
    run                   = psc_output_fields_c_run;
  }
} psc_output_fields_c_ops;
