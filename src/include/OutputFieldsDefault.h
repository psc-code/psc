
#pragma once

#include "../libpsc/psc_output_fields/fields_item_jeh.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "fields_item.hxx"

#include <mrc_io.hxx>

#include <memory>

// ======================================================================

using FieldsItem_E_cc = FieldsItemFields<ItemLoopPatches<Item_e_cc>>;
using FieldsItem_H_cc = FieldsItemFields<ItemLoopPatches<Item_h_cc>>;
using FieldsItem_J_cc = FieldsItemFields<ItemLoopPatches<Item_j_cc>>;

using _FieldsItem_n_1st_cc =
  _FieldsItemMoment<ItemMomentAddBnd<Moment_n_1st<MfieldsC>>>;
using _FieldsItem_v_1st_cc =
  _FieldsItemMoment<ItemMomentAddBnd<Moment_v_1st<MfieldsC>>>;
using _FieldsItem_p_1st_cc =
  _FieldsItemMoment<ItemMomentAddBnd<Moment_p_1st<MfieldsC>>>;
using _FieldsItem_T_1st_cc =
  _FieldsItemMoment<ItemMomentAddBnd<Moment_T_1st<MfieldsC>>>;

// ======================================================================
// OutputFieldsParams

struct OutputFieldsParams
{
  const char* data_dir = {"."};

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
// OutputFields

class OutputFields : public OutputFieldsParams
{
public:
  // ----------------------------------------------------------------------
  // ctor

  OutputFields(const Grid_t& grid, const OutputFieldsParams& prm)
    : OutputFieldsParams{prm},
      e_cc_{grid},
      h_cc_{grid},
      j_cc_{grid},
      n_cc_{grid},
      v_cc_{grid},
      T_cc_{grid}
  {
    pfield_next = pfield_first;
    tfield_next = tfield_first;

    tfds_.emplace_back(grid, e_cc_.n_comps(grid), grid.ibn);
    tfds_.emplace_back(grid, h_cc_.n_comps(grid), grid.ibn);
    tfds_.emplace_back(grid, j_cc_.n_comps(grid), grid.ibn);
    tfds_.emplace_back(grid, n_cc_.n_comps(grid), grid.ibn);
    tfds_.emplace_back(grid, v_cc_.n_comps(grid), grid.ibn);
    tfds_.emplace_back(grid, T_cc_.n_comps(grid), grid.ibn);

    naccum = 0;

    if (pfield_step > 0) {
      io_pfd_.reset(new MrcIo{"pfd", data_dir});
    }
    if (tfield_step) {
      io_tfd_.reset(new MrcIo{"tfd", data_dir});
    }
  }

  // ----------------------------------------------------------------------
  // operator()

  template <typename MfieldsState, typename Mparticles>
  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    const auto& grid = mflds._grid();

    static int pr;
    if (!pr) {
      pr = prof_register("output_c_field", 1., 0, 0);
    }
    prof_start(pr);

    auto timestep = grid.timestep();
    bool doaccum_tfield =
      tfield_step > 0 && (((timestep >= (tfield_next - tfield_length + 1)) &&
                           timestep % tfield_every == 0) ||
                          timestep == 0);

    if ((pfield_step > 0 && timestep >= pfield_next) ||
        (tfield_step > 0 && doaccum_tfield)) {
#if 0
	auto& mres = dynamic_cast<MfieldsC&>(item->mres());
	double min = 1e10, max = -1e10;
	for (int m = 0; m < mres.n_comps(); m++) {
	  for (int p = 0; p < mres.n_patches(); p++) {
	    mres.grid().Foreach_3d(0, 0, [&](int i, int j, int k) {
		min = fminf(min, mres[p](m, i,j,k));
		max = fmaxf(max, mres[p](m, i,j,k));
	      });
	  }
	  mprintf("name %s %d min %g max %g\n", item->name(), m, min, max);
	}
#endif
      e_cc_.run(grid, mflds, mprts);
      h_cc_.run(grid, mflds, mprts);
      j_cc_.run(grid, mflds, mprts);
      n_cc_.run(grid, mflds, mprts);
      v_cc_.run(grid, mflds, mprts);
      T_cc_.run(grid, mflds, mprts);
    }

    if (pfield_step > 0 && timestep >= pfield_next) {
      mpi_printf(MPI_COMM_WORLD, "***** Writing PFD output\n"); // FIXME
      pfield_next += pfield_step;

      io_pfd_->open(grid, rn, rx);
      write_pfd(e_cc_);
      write_pfd(h_cc_);
      write_pfd(j_cc_);
      write_pfd(n_cc_);
      write_pfd(v_cc_);
      write_pfd(T_cc_);
      io_pfd_->close();
    }

    if (tfield_step > 0) {
      if (doaccum_tfield) {
        // tfd += pfd
        int i = 0;
        tfds_[i++].axpy(1., e_cc_.mres());
        tfds_[i++].axpy(1., h_cc_.mres());
        tfds_[i++].axpy(1., j_cc_.mres());
        tfds_[i++].axpy(1., n_cc_.mres());
        tfds_[i++].axpy(1., v_cc_.mres());
        tfds_[i++].axpy(1., T_cc_.mres());
        naccum++;
      }
      if (timestep >= tfield_next) {
        mpi_printf(MPI_COMM_WORLD, "***** Writing TFD output\n"); // FIXME
        tfield_next += tfield_step;

        io_tfd_->open(grid, rn, rx);

        // convert accumulated values to correct temporal mean
        int i = 0;
        write_tfd(tfds_[i++], e_cc_);
        write_tfd(tfds_[i++], h_cc_);
        write_tfd(tfds_[i++], j_cc_);
        write_tfd(tfds_[i++], n_cc_);
        write_tfd(tfds_[i++], v_cc_);
        write_tfd(tfds_[i++], T_cc_);

        naccum = 0;
        io_tfd_->close();
      }
    }

    prof_stop(pr);
  };

private:
  template <typename Item>
  void write_pfd(Item& item)
  {
    item.mres().write_as_mrc_fld(io_pfd_->io_, item.name(), item.comp_names());
  }

  template <typename Item>
  void write_tfd(MfieldsC& tfd, Item& item)
  {
    tfd.scale(1. / naccum);
    tfd.write_as_mrc_fld(io_tfd_->io_, item.name(), item.comp_names());
    tfd.zero();
  }

public:
  int pfield_next, tfield_next;
  // storage for output
  unsigned int naccum;

private:
  FieldsItem_E_cc e_cc_;
  FieldsItem_H_cc h_cc_;
  FieldsItem_J_cc j_cc_;
  _FieldsItem_n_1st_cc n_cc_;
  _FieldsItem_n_1st_cc v_cc_;
  _FieldsItem_n_1st_cc T_cc_;
  // tfd -- FIXME?! always MfieldsC
  std::vector<MfieldsC> tfds_;
  std::unique_ptr<MrcIo> io_pfd_;
  std::unique_ptr<MrcIo> io_tfd_;
};

OutputFields defaultOutputFields(const Grid_t& grid,
                                 const OutputFieldsParams& params)
{
  return OutputFields{grid, params};
}
