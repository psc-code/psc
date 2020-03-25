
#pragma once

#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "fields_item.hxx"

#include <mrc_io.hxx>

#include <memory>

template <typename Mparticles>
using FieldsItem_Moments_1st_cc = Moments_1st<Mparticles>;

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
  using MfieldsFake = MfieldsC;
  using MparticlesFake = MparticlesDouble;

public:
  // ----------------------------------------------------------------------
  // ctor

  OutputFields(const Grid_t& grid, const OutputFieldsParams& prm)
    : OutputFieldsParams{prm},
      tfd_jeh_{grid, Item_jeh<MfieldsFake>::n_comps(), grid.ibn},
      tfd_moments_{grid, FieldsItem_Moments_1st_cc<MparticlesFake>::n_comps(grid), grid.ibn},
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

    if (!do_pfield && !doaccum_tfield) {
      prof_stop(pr);
      return;
    }

    Item_jeh<MfieldsState> pfd_jeh{mflds};
    FieldsItem_Moments_1st_cc<Mparticles> pfd_moments{mprts};

    if (do_pfield) {
      mpi_printf(grid.comm(), "***** Writing PFD output\n");
      pfield_next_ += pfield_step;

      io_pfd_->open(grid, rn, rx);
      _write_pfd(pfd_jeh);
      _write_pfd(pfd_moments);
      io_pfd_->close();
    }

    if (doaccum_tfield) {
      // tfd += pfd
      tfd_jeh_ += pfd_jeh;
      tfd_moments_ += pfd_moments;
      naccum_++;
    }
    if (do_tfield) {
      mpi_printf(grid.comm(), "***** Writing TFD output\n");
      tfield_next_ += tfield_step;

      io_tfd_->open(grid, rn, rx);
      _write_tfd(tfd_jeh_, pfd_jeh);
      _write_tfd(tfd_moments_, pfd_moments);
      io_tfd_->close();
      naccum_ = 0;
    }

    prof_stop(pr);
  };

private:
  template <typename EXP>
  void _write_pfd(EXP& pfd)
  {
    MrcIo::write_mflds(io_pfd_->io_, adaptMfields(pfd), pfd.grid(), pfd.name(),
                       pfd.comp_names());
  }

  template <typename EXP>
  void _write_tfd(MfieldsC& tfd, EXP& pfd)
  {
    // convert accumulated values to correct temporal mean
    tfd.scale(1. / naccum_);
    tfd.write_as_mrc_fld(io_tfd_->io_, pfd.name(), pfd.comp_names());
    tfd.zero();
  }

private:
  // tfd -- FIXME?! always MfieldsC
  MfieldsC tfd_jeh_;
  MfieldsC tfd_moments_;
  std::unique_ptr<MrcIo> io_pfd_;
  std::unique_ptr<MrcIo> io_tfd_;
  int pfield_next_, tfield_next_;
  int naccum_ = 0;
};
