
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
};

// ======================================================================
// OutputFieldsItem

template <typename Writer>
class OutputFieldsItem : public OutputFieldsItemParams
{
public:
  OutputFieldsItem(const Grid_t& grid, const OutputFieldsItemParams& prm,
                   int n_comps, Int3 ibn)
    : OutputFieldsItemParams{prm},
      pfield_next_{prm.pfield_first},
      tfield_next_{prm.tfield_first},
      tfd_{grid, n_comps, ibn}
  {}

  // private:
  int pfield_next_;
  int tfield_next_;
  Writer io_pfd_;
  Writer io_tfd_;
  // tfd -- FIXME?! always MfieldsC
  MfieldsC tfd_;
};

// ======================================================================
// OutputFieldsParams

struct OutputFieldsParams
{
  const char* data_dir = {"."};

  OutputFieldsItemParams fields;
  OutputFieldsItemParams moments;

  Int3 rn = {};
  Int3 rx = {1000000, 1000000, 100000};
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
      fields{grid, prm.fields, Item_jeh<MfieldsFake>::n_comps(), {}},
      moments{grid, prm.moments,
              FieldsItem_Moments_1st_cc<MparticlesFake>::n_comps(grid),
              grid.ibn},
      rn{prm.rn},
      rx{prm.rx}
  {
    if (fields.pfield_interval > 0) {
      fields.io_pfd_.open("pfd", data_dir);
    }
    if (moments.pfield_interval > 0) {
      moments.io_pfd_.open("pfd_moments", data_dir);
    }
    if (fields.tfield_interval > 0) {
      fields.io_tfd_.open("tfd", data_dir);
    }
    if (moments.tfield_interval > 0) {
      moments.io_tfd_.open("tfd_moments", data_dir);
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
        mpi_printf(grid.comm(), "***** Writing PFD output\n");
        fields.pfield_next_ += fields.pfield_interval;

        prof_start(pr_field_write);
        fields.io_pfd_.begin_step(grid);
        fields.io_pfd_.set_subset(grid, rn, rx);
        _write_pfd(fields.io_pfd_, pfd_jeh);
        fields.io_pfd_.end_step();
        prof_stop(pr_field_write);
      }

      if (doaccum_tfield) {
        // tfd += pfd
        prof_start(pr_field_acc);
        fields.tfd_ += pfd_jeh;
        prof_stop(pr_field_acc);
        naccum_++;
      }
      if (do_tfield) {
        mpi_printf(grid.comm(), "***** Writing TFD output\n");
        fields.tfield_next_ += fields.tfield_interval;

        prof_start(pr_field_write);
        fields.io_tfd_.begin_step(grid);
        fields.io_tfd_.set_subset(grid, rn, rx);
        _write_tfd(fields.io_tfd_, fields.tfd_, pfd_jeh, naccum_);
        fields.io_tfd_.end_step();
        naccum_ = 0;
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
        mpi_printf(grid.comm(), "***** Writing PFD moment output\n");
        moments.pfield_next_ += moments.pfield_interval;

        prof_start(pr_moment_write);
        moments.io_pfd_.begin_step(grid);
        moments.io_pfd_.set_subset(grid, rn, rx);
        _write_pfd(moments.io_pfd_, pfd_moments);
        moments.io_pfd_.end_step();
        prof_stop(pr_moment_write);
      }

      if (doaccum_tfield_moments) {
        // tfd += pfd
        prof_start(pr_moment_acc);
        moments.tfd_ += pfd_moments;
        prof_stop(pr_moment_acc);
        naccum_moments_++;
      }
      if (do_tfield_moments) {
        mpi_printf(grid.comm(), "***** Writing TFD moment output\n");
        moments.tfield_next_ += moments.tfield_interval;

        prof_start(pr_moment_write);
        moments.io_tfd_.begin_step(grid);
        moments.io_tfd_.set_subset(grid, rn, rx);
        _write_tfd(moments.io_tfd_, moments.tfd_, pfd_moments, naccum_moments_);
        moments.io_tfd_.end_step();
        prof_stop(pr_moment_write);
        naccum_moments_ = 0;
      }
      prof_stop(pr_moment);
    }

    prof_stop(pr);
  };

private:
  template <typename EXP>
  static void _write_pfd(Writer& io, EXP& pfd)
  {
    io.write(adapt(evalMfields(pfd)), pfd.grid(), pfd.name(), pfd.comp_names());
  }

  template <typename EXP>
  static void _write_tfd(Writer& io, MfieldsC& tfd, EXP& pfd, int naccum)
  {
    // convert accumulated values to correct temporal mean
    tfd.scale(1. / naccum);
    io.write(adapt(tfd), tfd.grid(), pfd.name(), pfd.comp_names());
    tfd.zero();
  }

public:
  const char* data_dir;
  OutputFieldsItem<Writer> fields;
  OutputFieldsItem<Writer> moments;
  Int3 rn = {};
  Int3 rx = {1000000, 1000000, 100000};

private:
  int naccum_ = 0;
  int naccum_moments_ = 0;
  bool first_time = true;
};

#ifdef xPSC_HAVE_ADIOS2

#include "writer_adios2.hxx"
using OutputFields = OutputFieldsDefault<WriterADIOS2>;

#else

#include "writer_mrc.hxx"
using OutputFields = OutputFieldsDefault<WriterMRC>;

#endif
