
#pragma once

#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "fields_item.hxx"
#include "psc_particles_double.h"

#include <memory>

template <typename Mparticles>
using FieldsItem_Moments_1st_cc = Moments_1st<Mparticles>;

// ======================================================================
// OutputFieldsItemParams

struct OutputFieldsItemParams
{
  int pfield_interval = 0;
  int pfield_first = 0;
  int tfield_interval = 0;
  int tfield_first = 0;
  int tfield_average_length = 1000000;
  int tfield_average_every = 1;
  Int3 rn = {};
  Int3 rx = {1000000, 1000000, 100000};
};

// ======================================================================
// OutputFieldsItem

template <typename Writer>
class OutputFieldsItem : public OutputFieldsItemParams
{
public:
  OutputFieldsItem(const Grid_t& grid, const OutputFieldsItemParams& prm,
                   int n_comps, Int3 ibn, const char* data_dir, std::string sfx)
    : OutputFieldsItemParams{prm},
      pfield_next_{prm.pfield_first},
      tfield_next_{prm.tfield_first},
      tfd_{grid, n_comps, ibn}
  {
    if (pfield_interval > 0) {
      io_pfd_.open("pfd" + sfx, data_dir);
    }
    if (tfield_interval > 0) {
      io_tfd_.open("tfd" + sfx, data_dir);
    }
  }

  template <typename EXP>
  void write_pfd(EXP& pfd)
  {
    mpi_printf(pfd.grid().comm(), "***** Writing PFD output\n");
    pfield_next_ += pfield_interval;
    io_pfd_.begin_step(pfd.grid());
    io_pfd_.set_subset(pfd.grid(), rn, rx);
    io_pfd_.write(adapt(evalMfields(pfd)), pfd.grid(), pfd.name(),
                  pfd.comp_names());
    io_pfd_.end_step();
  }

  template <typename EXP>
  void write_tfd(MfieldsC& tfd, EXP& pfd)
  {
    mpi_printf(tfd.grid().comm(), "***** Writing TFD output\n");
    tfield_next_ += tfield_interval;

    // convert accumulated values to correct temporal mean
    tfd.scale(1. / naccum_);

    io_tfd_.begin_step(tfd.grid());
    io_tfd_.set_subset(tfd.grid(), rn, rx);
    io_tfd_.write(adapt(tfd), tfd.grid(), pfd.name(), pfd.comp_names());
    io_tfd_.end_step();
    naccum_ = 0;

    tfd.zero();
  }

  // private:
  int pfield_next_;
  int tfield_next_;
  Writer io_pfd_;
  Writer io_tfd_;
  // tfd -- FIXME?! always MfieldsC
  MfieldsC tfd_;
  int naccum_ = 0;
};

// ======================================================================
// OutputFieldsParams

struct OutputFieldsParams
{
  const char* data_dir = {"."};

  OutputFieldsItemParams fields;
  OutputFieldsItemParams moments;
};

// ======================================================================
// OutputFieldsDefault

template <typename Writer>
class OutputFieldsDefault
{
  using MfieldsFake = MfieldsC;
  using MparticlesFake = MparticlesDouble;

public:
  // ----------------------------------------------------------------------
  // ctor

  OutputFieldsDefault(const Grid_t& grid, const OutputFieldsParams& prm)
    : data_dir{prm.data_dir},
      fields{grid, prm.fields, Item_jeh<MfieldsFake>::n_comps(),
             {},   data_dir,   ""},
      moments{grid,
              prm.moments,
              FieldsItem_Moments_1st_cc<MparticlesFake>::n_comps(grid),
              grid.ibn,
              data_dir,
              "_moments"}
  {}

  // ----------------------------------------------------------------------
  // operator()

  template <typename MfieldsState, typename Mparticles>
  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    const auto& grid = mflds._grid();

    static int pr;
    if (!pr) {
      pr = prof_register("outf", 1., 0, 0);
    }
#if 1
    static int pr_field, pr_moment, pr_field_calc, pr_moment_calc,
      pr_field_write, pr_moment_write, pr_field_acc, pr_moment_acc;
    if (!pr_field) {
      pr_field = prof_register("outf_field", 1., 0, 0);
      pr_moment = prof_register("outf_moment", 1., 0, 0);
      pr_field_calc = prof_register("outf_field_calc", 1., 0, 0);
      pr_moment_calc = prof_register("outf_moment_calc", 1., 0, 0);
      pr_field_write = prof_register("outf_field_write", 1., 0, 0);
      pr_moment_write = prof_register("outf_moment_write", 1., 0, 0);
      pr_field_acc = prof_register("outf_field_acc", 1., 0, 0);
      pr_moment_acc = prof_register("outf_moment_acc", 1., 0, 0);
    }
#endif

    auto timestep = grid.timestep();
    if (first_time) {
      first_time = false;
      if (timestep != 0) {
        fields.pfield_next_ = timestep + fields.pfield_interval;
        fields.tfield_next_ = timestep + fields.tfield_interval;
        moments.pfield_next_ = timestep + moments.pfield_interval;
        moments.tfield_next_ = timestep + moments.tfield_interval;
        return;
      }
    }

    prof_start(pr);

    bool do_pfield =
      fields.pfield_interval > 0 && timestep >= fields.pfield_next_;
    bool do_tfield =
      fields.tfield_interval > 0 && timestep >= fields.tfield_next_;
    bool doaccum_tfield = fields.tfield_interval > 0 &&
                          (((timestep >= (fields.tfield_next_ -
                                          fields.tfield_average_length + 1)) &&
                            timestep % fields.tfield_average_every == 0) ||
                           timestep == 0);

    if (do_pfield || doaccum_tfield) {
      prof_start(pr_field);
      prof_start(pr_field_calc);
      Item_jeh<MfieldsState> pfd_jeh{mflds};
      prof_stop(pr_field_calc);

      if (do_pfield) {
        prof_start(pr_field_write);
        fields.write_pfd(pfd_jeh);
        prof_stop(pr_field_write);
      }

      if (doaccum_tfield) {
        // tfd += pfd
        prof_start(pr_field_acc);
        fields.tfd_ += pfd_jeh;
        prof_stop(pr_field_acc);
        fields.naccum_++;
      }
      if (do_tfield) {
        prof_start(pr_field_write);
        fields.write_tfd(fields.tfd_, pfd_jeh);
        prof_stop(pr_field_write);
      }
      prof_stop(pr_field);
    }

    bool do_pfield_moments =
      moments.pfield_interval > 0 && timestep >= moments.pfield_next_;
    bool do_tfield_moments =
      moments.tfield_interval > 0 && timestep >= moments.tfield_next_;
    bool doaccum_tfield_moments =
      moments.tfield_interval > 0 &&
      (((timestep >=
         (moments.tfield_next_ - moments.tfield_average_length + 1)) &&
        timestep % moments.tfield_average_every == 0) ||
       timestep == 0);

    if (do_pfield_moments || doaccum_tfield_moments) {
      prof_start(pr_moment);
      prof_start(pr_moment_calc);
      FieldsItem_Moments_1st_cc<Mparticles> pfd_moments{mprts};
      prof_stop(pr_moment_calc);

      if (do_pfield_moments) {
        prof_start(pr_moment_write);
        moments.write_pfd(pfd_moments);
        prof_stop(pr_moment_write);
      }

      if (doaccum_tfield_moments) {
        // tfd += pfd
        prof_start(pr_moment_acc);
        moments.tfd_ += pfd_moments;
        prof_stop(pr_moment_acc);
        moments.naccum_++;
      }
      if (do_tfield_moments) {
        prof_start(pr_moment_write);
        moments.write_tfd(moments.tfd_, pfd_moments);
        prof_stop(pr_moment_write);
      }
      prof_stop(pr_moment);
    }

    prof_stop(pr);
  };

public:
  const char* data_dir;
  OutputFieldsItem<Writer> fields;
  OutputFieldsItem<Writer> moments;

private:
  bool first_time = true;
};

#ifdef xPSC_HAVE_ADIOS2

#include "writer_adios2.hxx"
using OutputFields = OutputFieldsDefault<WriterADIOS2>;

#else

#include "writer_mrc.hxx"
using OutputFields = OutputFieldsDefault<WriterMRC>;

#endif
