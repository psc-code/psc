
#pragma once

#include "../libpsc/psc_output_fields/fields_item_jeh.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "fields_item.hxx"

#include <mrc_io.hxx>

#include <memory>

// ======================================================================

using FieldsItem_jeh = FieldsItemFields<ItemLoopPatches<Item_jeh>>;
using FieldsItem_E_cc = FieldsItemFields<ItemLoopPatches<Item_e_cc>>;
using FieldsItem_H_cc = FieldsItemFields<ItemLoopPatches<Item_h_cc>>;
using FieldsItem_J_cc = FieldsItemFields<ItemLoopPatches<Item_j_cc>>;

using FieldsItem_n_1st_cc =
  FieldsItemMoment<ItemMomentAddBnd<Moment_n_1st<MfieldsC>>>;
using FieldsItem_v_1st_cc =
  FieldsItemMoment<ItemMomentAddBnd<Moment_v_1st<MfieldsC>>>;
using FieldsItem_p_1st_cc =
  FieldsItemMoment<ItemMomentAddBnd<Moment_p_1st<MfieldsC>>>;
using FieldsItem_T_1st_cc =
  FieldsItemMoment<ItemMomentAddBnd<Moment_T_1st<MfieldsC>>>;
using FieldsItem_Moments_1st_cc =
  FieldsItemMoment<ItemMomentAddBnd<Moments_1st<MfieldsC>>>;

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
      jeh_{grid},
      moments_{grid},
      tfd_jeh_{grid, jeh_.n_comps(grid), grid.ibn},
      tfd_moments_{grid, moments_.n_comps(grid), grid.ibn},
      pfield_next_{pfield_first},
      tfield_next_{tfield_first}
  {
    if (pfield_step > 0) {
      io_pfd_.reset(new MrcIo{"pfd", data_dir});
    }
    if (tfield_step > 0) {
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
    bool do_pfield = pfield_step > 0 && timestep >= pfield_next_;
    bool do_tfield = tfield_step > 0 && timestep >= tfield_next_;
    bool doaccum_tfield =
      tfield_step > 0 && (((timestep >= (tfield_next_ - tfield_length + 1)) &&
                           timestep % tfield_every == 0) ||
                          timestep == 0);

    if ((do_pfield || doaccum_tfield)) {
      jeh_.run(grid, mflds, mprts);
      moments_.run(grid, mflds, mprts);
    }

    if (do_pfield) {
      mpi_printf(grid.comm(), "***** Writing PFD output\n");
      pfield_next_ += pfield_step;

      io_pfd_->open(grid, rn, rx);
      write_pfd(jeh_);
      write_pfd(moments_);
      io_pfd_->close();
    }

    if (doaccum_tfield) {
      // tfd += pfd
      tfd_jeh_.axpy(1., jeh_.mres());
      tfd_moments_.axpy(1., moments_.mres());
      naccum_++;
    }
    if (do_tfield) {
      mpi_printf(grid.comm(), "***** Writing TFD output\n");
      tfield_next_ += tfield_step;
      
      io_tfd_->open(grid, rn, rx);
      write_tfd(tfd_jeh_, jeh_);
      write_tfd(tfd_moments_, moments_);
      io_tfd_->close();
      naccum_ = 0;
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
    // convert accumulated values to correct temporal mean
    tfd.scale(1. / naccum_);
    tfd.write_as_mrc_fld(io_tfd_->io_, item.name(), item.comp_names());
    tfd.zero();
  }

private:
  FieldsItem_jeh jeh_;
  FieldsItem_Moments_1st_cc moments_;
  // tfd -- FIXME?! always MfieldsC
  MfieldsC tfd_jeh_;
  MfieldsC tfd_moments_;
  std::unique_ptr<MrcIo> io_pfd_;
  std::unique_ptr<MrcIo> io_tfd_;
  int pfield_next_, tfield_next_;
  int naccum_ = 0;
};

