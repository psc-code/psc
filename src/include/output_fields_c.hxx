
#pragma once

#include "fields_item.hxx"

#include <mrc_io.hxx>

// ======================================================================
// OutputFieldsCParams

struct OutputFieldsCParams
{
  const char *data_dir = {"."};
  const char *output_fields = {"j,e,h"};

  int pfield_step = 0;
  int pfield_first = 0;

  int tfield_step = 0;
  int tfield_first = 0;
  int tfield_length = 1000000;
  int tfield_every = 1;

  Int3 rn = {};
  Int3 rx = {1000000, 1000000, 100000};
};

// ======================================================================
// OutputFieldsC

struct OutputFieldsC : public OutputFieldsCParams
{
  struct Item
  {
    Item(PscFieldsItemBase item, const std::string& name,
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

  // ----------------------------------------------------------------------
  // ctor

  OutputFieldsC(const Grid_t& grid, const OutputFieldsCParams& prm)
    : OutputFieldsCParams{prm}
  {
    pfield_next = pfield_first;
    tfield_next = tfield_first;

    if (output_fields) {
      // setup pfd according to output_fields as given
      // (potentially) on the command line
      // parse comma separated list of fields
      char *s_orig = strdup(output_fields), *p, *s = s_orig;
      while ((p = strsep(&s, ", "))) {
	struct psc_output_fields_item *item =
	  psc_output_fields_item_create(grid.comm());
	psc_output_fields_item_set_type(item, p);
	psc_output_fields_item_setup(item);
	
	// pfd
	std::vector<std::string> comp_names = PscFieldsItemBase{item}->comp_names();
	MfieldsBase& mflds_pfd = PscFieldsItemBase{item}->mres();
	
	// tfd -- FIXME?! always MfieldsC
	MfieldsBase& mflds_tfd = *new MfieldsC{grid, mflds_pfd.n_comps(), grid.ibn};
	items.emplace_back(PscFieldsItemBase{item}, p, comp_names, mflds_pfd, mflds_tfd);
      }
      free(s_orig);
    }
    
    naccum = 0;
    
    if (pfield_step > 0) {
      io_pfd_.reset(new MrcIo{"pfd", data_dir});
    }
    if (tfield_step) {
      io_tfd_.reset(new MrcIo{"tfd", data_dir});
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
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsStateBase& mflds, MparticlesBase& mprts)
  {
    const auto& grid = mflds.grid();
    
    static int pr;
    if (!pr) {
      pr = prof_register("output_c_field", 1., 0, 0);
    }
    prof_start(pr);

    auto timestep = grid.timestep();
    bool doaccum_tfield = tfield_step > 0 && 
      (((timestep >= (tfield_next - tfield_length + 1)) &&
	timestep % tfield_every == 0) ||
       timestep == 0);
    
    if ((pfield_step > 0 && timestep >= pfield_next) ||
	(tfield_step > 0 && doaccum_tfield)) {
      for (auto item : items) {
	item.item(mflds, mprts);
      }
    }
    
    if (pfield_step > 0 && timestep >= pfield_next) {
      mpi_printf(MPI_COMM_WORLD, "***** Writing PFD output\n"); // FIXME
      pfield_next += pfield_step;
      
      io_pfd_->open(grid, rn, rx);
      for (auto& item : items) {
	item.pfd.write_as_mrc_fld(io_pfd_->io_, item.name, item.comp_names);
      }
      io_pfd_->close();
    }
    
    if (tfield_step > 0) {
      if (doaccum_tfield) {
	// tfd += pfd
	for (auto& item : items) {
	  item.tfd.axpy(1., item.pfd);
	}
	naccum++;
      }
      if (timestep >= tfield_next) {
	mpi_printf(MPI_COMM_WORLD, "***** Writing TFD output\n"); // FIXME
	tfield_next += tfield_step;
	
	io_tfd_->open(grid, rn, rx);
	
	// convert accumulated values to correct temporal mean
	for (auto& item : items) {
	  item.tfd.scale(1. / naccum);
	  item.tfd.write_as_mrc_fld(io_tfd_->io_, item.name, item.comp_names);
	  item.tfd.zero();
	}
	naccum = 0;
	io_tfd_->close();
      }
    }

    prof_stop(pr);
  };

public:
  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;
  std::vector<Item> items;
private:
  std::unique_ptr<MrcIo> io_pfd_;
  std::unique_ptr<MrcIo> io_tfd_;
};

