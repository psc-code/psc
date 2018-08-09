
#include "psc_output_fields_private.h"
#include "psc_output_fields_item.h"
#include "psc_fields_as_c.h"
#include "bnd.hxx"
#include "fields_item.hxx"
#include "output_fields.hxx"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <string.h>


#include "fields3d.hxx"
#include "fields_item.hxx"

// ----------------------------------------------------------------------

struct OutputFieldsItem
{
  OutputFieldsItem(PscFieldsItemBase item, const std::string& name,
		   std::vector<std::string>& comp_names, MfieldsBase& pfd,
		   MfieldsBase& tfd)
    : item(item), name(name), comp_names(comp_names), pfd(pfd), tfd(tfd)
  {}

  PscFieldsItemBase item;
  MfieldsBase& pfd;
  MfieldsBase& tfd;
  std::string name;
  std::vector<std::string> comp_names;
};

// ======================================================================
// OutputFieldsCParams

struct OutputFieldsCParams
{
  char *data_dir;
  char *output_fields;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_first, tfield_first;
  int pfield_step, tfield_step;
  int tfield_length;
  int tfield_every;
  int rn[3];
  int rx[3];
};

// ======================================================================
// OutputFieldsC

struct OutputFieldsC : public OutputFieldsCParams
{
  // ----------------------------------------------------------------------
  // ctor

  OutputFieldsC(MPI_Comm comm, const OutputFieldsCParams& prm)
    : OutputFieldsCParams{prm}
  {
    pfield_next = pfield_first;
    tfield_next = tfield_first;

    struct psc *psc = ppsc;
    
    // setup pfd according to output_fields as given
    // (potentially) on the command line
    // parse comma separated list of fields
    char *s_orig = strdup(output_fields), *p, *s = s_orig;
    while ((p = strsep(&s, ", "))) {
      struct psc_output_fields_item *item =
	psc_output_fields_item_create(comm);
      psc_output_fields_item_set_type(item, p);
      psc_output_fields_item_setup(item);
      
      // pfd
      std::vector<std::string> comp_names = PscFieldsItemBase{item}->comp_names();
      MfieldsBase& mflds_pfd = PscFieldsItemBase{item}->mres();
      
      // tfd -- FIXME?! always MfieldsC
      MfieldsBase& mflds_tfd = *new MfieldsC{psc->grid(), mflds_pfd.n_comps(), psc->ibn};
      items.emplace_back(PscFieldsItemBase{item}, p, comp_names, mflds_pfd, mflds_tfd);
    }
    free(s_orig);
    
    naccum = 0;
    
    if (dowrite_pfield) {
      io_pfd_ = create_mrc_io("pfd");
    }
    if (dowrite_tfield) {
      io_tfd_ = create_mrc_io("tfd");
    }
  }

  // ----------------------------------------------------------------------
  // dtor

  ~OutputFieldsC()
  {
    for (auto& item : items) {
      psc_output_fields_item_destroy(item.item.item());
      delete &item.tfd;
    }
    
    mrc_io_destroy(io_pfd_);
    mrc_io_destroy(io_tfd_);
  }

  // ----------------------------------------------------------------------
  // run

  void run(MfieldsBase& mflds, MparticlesBase& mprts)
  {
    auto psc = ppsc;
    
    static int pr;
    if (!pr) {
      pr = prof_register("output_c_field", 1., 0, 0);
    }
    prof_start(pr);

    bool doaccum_tfield = dowrite_tfield && 
      (((psc->timestep >= (tfield_next - tfield_length + 1)) &&
	psc->timestep % tfield_every == 0) ||
       psc->timestep == 0);
    
    if ((dowrite_pfield && psc->timestep >= pfield_next) ||
	(dowrite_tfield && doaccum_tfield)) {
      for (auto item : items) {
	item.item(mflds, mprts);
      }
    }
    
    if (dowrite_pfield && psc->timestep >= pfield_next) {
      mpi_printf(MPI_COMM_WORLD, "***** Writing PFD output\n"); // FIXME
      pfield_next += pfield_step;
      
      open_mrc_io(io_pfd_);
      for (auto& item : items) {
	item.pfd.write_as_mrc_fld(io_pfd_, item.name, item.comp_names);
      }
      mrc_io_close(io_pfd_);
    }
    
    if (dowrite_tfield) {
      if (doaccum_tfield) {
	// tfd += pfd
	for (auto& item : items) {
	  item.tfd.axpy(1., item.pfd);
	}
	naccum++;
      }
      if (psc->timestep >= tfield_next) {
	mpi_printf(MPI_COMM_WORLD, "***** Writing TFD output\n"); // FIXME
	tfield_next += tfield_step;
	
	open_mrc_io(io_tfd_);
	
	// convert accumulated values to correct temporal mean
	for (auto& item : items) {
	  item.tfd.scale(1. / naccum);
	  item.tfd.write_as_mrc_fld(io_tfd_, item.name, item.comp_names);
	  item.tfd.zero();
	}
	naccum = 0;
	mrc_io_close(io_tfd_);
      }
    }

    prof_stop(pr);
  };

private:
  mrc_io* create_mrc_io(const char* pfx)
  {
    mrc_io* io = mrc_io_create(MPI_COMM_WORLD);
    mrc_io_set_param_string(io, "basename", pfx);
    mrc_io_set_param_string(io, "outdir", data_dir);
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
    mrc_io_view(io);
    return io;
  }

  void open_mrc_io(mrc_io *io)
  {
    int gdims[3];
    mrc_domain_get_global_dims(ppsc->mrc_domain_, gdims);
    int slab_off[3], slab_dims[3];
    for (int d = 0; d < 3; d++) {
      if (rx[d] > gdims[d])
	rx[d] = gdims[d];
      
      slab_off[d] = rn[d];
      slab_dims[d] = rx[d] - rn[d];
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

public:
  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  std::vector<OutputFieldsItem> items;
private:
  mrc_io *io_pfd_;
  mrc_io *io_tfd_;

};

using PscOutputFields_t = PscOutputFields<OutputFieldsC>;

// ----------------------------------------------------------------------
// OutputFieldsC_destroy

static void
OutputFieldsC_destroy(struct psc_output_fields *out)
{
  PscOutputFields_t outf{out};

  outf->~OutputFieldsC();
}

// ----------------------------------------------------------------------
// OutputFieldsC_setup

static void
OutputFieldsC_setup(struct psc_output_fields *out)
{
  PscOutputFields_t outf{out};

  OutputFieldsCParams prm = *outf.sub();
  new(outf.sub()) OutputFieldsC{psc_output_fields_comm(out), prm};
}

// ----------------------------------------------------------------------
// OutputFieldsC_run

static void
OutputFieldsC_run(struct psc_output_fields *out,
			MfieldsBase& mflds, MparticlesBase& mprts)
{
  PscOutputFields_t outf{out};

  outf->run(mflds, mprts);
}

// ======================================================================
// psc_output_fields: subclass "c"

#define VAR(x) (void *)offsetof(struct OutputFieldsC, x)

// FIXME pfield_out_[xyz]_{min,max} aren't for pfield only, better init to 0,
// use INT3

static struct param OutputFieldsC_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "output_fields"      , VAR(output_fields)        , PARAM_STRING("j,e,h")   },
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
    size                  = sizeof(struct OutputFieldsC);
    param_descr           = OutputFieldsC_descr;
    setup                 = OutputFieldsC_setup;
    destroy               = OutputFieldsC_destroy;
    run                   = OutputFieldsC_run;
  }
} psc_output_fields_c_ops;
