
#pragma once

#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#include "fields_item.hxx"
#include "psc_particles_double.h"

#include <memory>

template <typename Mparticles>
using FieldsItem_Moments_1st_cc = Moments_1st<Mparticles>;

// ======================================================================
// OutputFieldsParams

struct OutputFieldsParams
{
  const char* data_dir = {"."};

  int pfield_interval = 0;
  int pfield_first = 0;

  int pfield_moments_interval = -1;
  int pfield_moments_first = -1;

  int tfield_interval = 0;
  int tfield_first = 0;
  int tfield_average_length = 1000000;
  int tfield_average_every = 1;

  int tfield_moments_interval = -1;
  int tfield_moments_first = -1;
  int tfield_moments_average_length = -1;
  int tfield_moments_average_every = -1;

  Int3 rn = {};
  Int3 rx = {1000000, 1000000, 100000};
};

// ======================================================================
// OutputFieldsDefault

template <typename Writer>
class OutputFieldsDefault : public OutputFieldsParams
{
  using MfieldsFake = MfieldsC;
  using MparticlesFake = MparticlesDouble;

public:
  // ----------------------------------------------------------------------
  // ctor

  OutputFieldsDefault(const Grid_t& grid, const OutputFieldsParams& prm)
    : OutputFieldsParams{prm},
      tfd_jeh_{grid, Item_jeh<MfieldsFake>::n_comps(), {}},
      tfd_moments_{grid,
                   FieldsItem_Moments_1st_cc<MparticlesFake>::n_comps(grid),
                   grid.ibn},
      pfield_next_{pfield_first},
      tfield_next_{tfield_first}
  {
    if (pfield_moments_interval < 0) {
      pfield_moments_interval = pfield_interval;
    }
    if (pfield_moments_first < 0) {
      pfield_moments_first = pfield_first;
    }

    if (tfield_moments_interval < 0) {
      tfield_moments_interval = tfield_interval;
    }
    if (tfield_moments_first < 0) {
      tfield_moments_first = tfield_first;
    }
    if (tfield_moments_average_length < 0) {
      tfield_moments_average_length = tfield_average_length;
    }
    if (tfield_moments_average_every < 0) {
      tfield_moments_average_every = tfield_average_every;
    }

    pfield_moments_next_ = pfield_moments_first;
    tfield_moments_next_ = tfield_moments_first;

    if (pfield_interval > 0) {
      io_pfd_.open("pfd", data_dir);
    }
    if (pfield_moments_interval > 0) {
      io_pfd_moments_.open("pfd_moments", data_dir);
    }
    if (tfield_interval > 0) {
      io_tfd_.open("tfd", data_dir);
    }
    if (tfield_moments_interval > 0) {
      io_tfd_moments_.open("tfd_moments", data_dir);
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
        pfield_next_ = timestep + pfield_interval;
        tfield_next_ = timestep + tfield_interval;
        pfield_moments_next_ = timestep + pfield_moments_interval;
        tfield_moments_next_ = timestep + tfield_moments_interval;
        return;
      }
    }

    prof_start(pr);

    bool do_pfield = pfield_interval > 0 && timestep >= pfield_next_;
    bool do_tfield = tfield_interval > 0 && timestep >= tfield_next_;
    bool doaccum_tfield =
      tfield_interval > 0 &&
      (((timestep >= (tfield_next_ - tfield_average_length + 1)) &&
        timestep % tfield_average_every == 0) ||
       timestep == 0);

    if (do_pfield || doaccum_tfield) {
      prof_start(pr_field);
      prof_start(pr_field_calc);
      Item_jeh<MfieldsState> pfd_jeh{mflds};
      prof_stop(pr_field_calc);

      if (do_pfield) {
        mpi_printf(grid.comm(), "***** Writing PFD output\n");
        pfield_next_ += pfield_interval;

        prof_start(pr_field_write);
        io_pfd_.begin_step(grid);
        io_pfd_.set_subset(grid, rn, rx);
        _write_pfd(io_pfd_, pfd_jeh);
        io_pfd_.end_step();
        prof_stop(pr_field_write);
      }

      if (doaccum_tfield) {
        // tfd += pfd
        prof_start(pr_field_acc);
        tfd_jeh_ += pfd_jeh;
        prof_stop(pr_field_acc);
        naccum_++;
      }
      if (do_tfield) {
        mpi_printf(grid.comm(), "***** Writing TFD output\n");
        tfield_next_ += tfield_interval;

        prof_start(pr_field_write);
        io_tfd_.begin_step(grid);
        io_tfd_.set_subset(grid, rn, rx);
        _write_tfd(io_tfd_, tfd_jeh_, pfd_jeh, naccum_);
        io_tfd_.end_step();
        naccum_ = 0;
        prof_stop(pr_field_write);
      }
      prof_stop(pr_field);
    }

    bool do_pfield_moments =
      pfield_moments_interval > 0 && timestep >= pfield_moments_next_;
    bool do_tfield_moments =
      tfield_moments_interval > 0 && timestep >= tfield_moments_next_;
    bool doaccum_tfield_moments =
      tfield_moments_interval > 0 &&
      (((timestep >=
         (tfield_moments_next_ - tfield_moments_average_length + 1)) &&
        timestep % tfield_moments_average_every == 0) ||
       timestep == 0);

    if (do_pfield_moments || doaccum_tfield_moments) {
      prof_start(pr_moment);
      prof_start(pr_moment_calc);
      FieldsItem_Moments_1st_cc<Mparticles> pfd_moments{mprts};
      prof_stop(pr_moment_calc);

      if (do_pfield_moments) {
        mpi_printf(grid.comm(), "***** Writing PFD moment output\n");
        pfield_moments_next_ += pfield_moments_interval;

        prof_start(pr_moment_write);
        io_pfd_moments_.begin_step(grid);
        io_pfd_moments_.set_subset(grid, rn, rx);
        _write_pfd(io_pfd_moments_, pfd_moments);
        io_pfd_moments_.end_step();
        prof_stop(pr_moment_write);
      }

      if (doaccum_tfield_moments) {
        // tfd += pfd
        prof_start(pr_moment_acc);
        tfd_moments_ += pfd_moments;
        prof_stop(pr_moment_acc);
        naccum_moments_++;
      }
      if (do_tfield_moments) {
        mpi_printf(grid.comm(), "***** Writing TFD moment output\n");
        tfield_moments_next_ += tfield_moments_interval;

        prof_start(pr_moment_write);
        io_tfd_moments_.begin_step(grid);
        io_tfd_moments_.set_subset(grid, rn, rx);
        _write_tfd(io_tfd_moments_, tfd_moments_, pfd_moments, naccum_moments_);
        io_tfd_moments_.end_step();
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
    io.write(evalMfields(pfd), pfd.grid(), pfd.name(), pfd.comp_names());
  }

  template <typename EXP>
  static void _write_tfd(Writer& io, MfieldsC& tfd, EXP& pfd, int naccum)
  {
    // convert accumulated values to correct temporal mean
    tfd.scale(1. / naccum);
    io.write(tfd, tfd.grid(), pfd.name(), pfd.comp_names());
    tfd.zero();
  }

private:
  // tfd -- FIXME?! always MfieldsC
  MfieldsC tfd_jeh_;
  MfieldsC tfd_moments_;
  Writer io_pfd_;
  Writer io_pfd_moments_;
  Writer io_tfd_;
  Writer io_tfd_moments_;
  int pfield_next_;
  int pfield_moments_next_;
  int tfield_next_;
  int tfield_moments_next_;
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
