
#pragma once

#include "fields_item.hxx"
#include "../libpsc/psc_output_fields/fields_item_jeh.hxx"

#include <mrc_io.hxx>

#include <memory>

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
	auto item = FieldsItemFactory::create(p, grid);
	
	// tfd -- FIXME?! always MfieldsC
	tfds_.emplace_back(new MfieldsC{grid, item->n_comps(grid), grid.ibn});
	items_.emplace_back(std::move(item));
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
    for (auto tfd : tfds_) {
      delete tfd;
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
      for (auto& item : items_) {
	item->run(mflds, mprts);
      }
    }
    
    if (pfield_step > 0 && timestep >= pfield_next) {
      mpi_printf(MPI_COMM_WORLD, "***** Writing PFD output\n"); // FIXME
      pfield_next += pfield_step;
      
      io_pfd_->open(grid, rn, rx);
      for (auto& item : items_) {
	item->mres().write_as_mrc_fld(io_pfd_->io_, item->_name(), item->comp_names());
      }
      io_pfd_->close();
    }
    
    if (tfield_step > 0) {
      if (doaccum_tfield) {
	// tfd += pfd
	for (int i = 0; i < items_.size(); i++) {
	  tfds_[i]->axpy(1., items_[i]->mres());
	}
	naccum++;
      }
      if (timestep >= tfield_next) {
	mpi_printf(MPI_COMM_WORLD, "***** Writing TFD output\n"); // FIXME
	tfield_next += tfield_step;
	
	io_tfd_->open(grid, rn, rx);
	
	// convert accumulated values to correct temporal mean
	for (int i = 0; i < items_.size(); i++) {
	  tfds_[i]->scale(1. / naccum);
	  tfds_[i]->write_as_mrc_fld(io_tfd_->io_, items_[i]->_name(), items_[i]->comp_names());
	  tfds_[i]->zero();
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
private:
  std::vector<std::unique_ptr<FieldsItemBase>> items_;
  std::vector<MfieldsBase*> tfds_;
  std::unique_ptr<MrcIo> io_pfd_;
  std::unique_ptr<MrcIo> io_tfd_;
};

