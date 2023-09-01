
#pragma once

#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "../libpsc/psc_output_fields/fields_item_moments_1st.hxx"
#ifdef USE_CUDA
#include "../libpsc/cuda/fields_item_moments_1st_cuda.hxx"
#endif
#include "fields_item.hxx"
#include "psc_particles_double.h"

#ifdef PSC_HAVE_ADIOS2

#include "writer_adios2.hxx"
using WriterDefault = WriterADIOS2;

#else

#include "writer_mrc.hxx"
using WriterDefault = WriterMRC;

#endif

#include <memory>

namespace detail
{
template <typename S, typename D, typename Enable = void>
struct moment_selector
{
  using type = Moments_1st<S, D>;
};

#ifdef USE_CUDA
template <typename S, typename D>
struct moment_selector<S, D,
                       std::enable_if_t<std::is_same<typename S::space_type,
                                                     gt::space::device>::value>>
{
  using type = Moments_1st_cuda<D>;
};
#endif
} // namespace detail

template <typename S, typename D>
using Item_Moments = typename detail::moment_selector<S, D>::type;

template <typename E>
struct mfields_gt
{
  E gt;
  std::string name;
  std::vector<std::string> comp_names;
};

template <typename E>
auto make_mfields_gt(E&& gt, const std::string& name,
                     const std::vector<std::string>& comp_names)
{
  return mfields_gt<E>{std::forward<E>(gt), name, comp_names};
}

// ======================================================================
// OutputFieldsItemParams

struct OutputFieldsItemParams
{
  std::string data_dir = ".";
  // distance between timesteps at which pfields are output (0 = disable)
  int pfield_interval = 0;
  // first timestep at which pfields are output
  int pfield_first = 0;
  // distance between timesteps at which tfields are output (0 = disable)
  int tfield_interval = 0;
  // first timestep at which tfields are output
  int tfield_first = 0;
  // max range of timesteps over which to average (capped at `tfield_interval`)
  int tfield_average_length = 1000000;
  // difference between timesteps used for average
  int tfield_average_every = 1;
  Int3 rn = {};
  Int3 rx = {10000000, 10000000, 10000000};

  bool tfield_enabled() { return tfield_interval > 0; }
};

// ======================================================================
// OutputFieldsItem

template <typename Mfields, typename Writer>
class OutputFieldsItem : public OutputFieldsItemParams
{
public:
  OutputFieldsItem(const Grid_t& grid, const OutputFieldsItemParams& prm,
                   std::string sfx)
    : OutputFieldsItemParams{prm},
      pfield_next_{prm.pfield_first},
      tfield_next_{prm.tfield_first}
  {
    if (pfield_interval > 0) {
      io_pfd_.open("pfd" + sfx, prm.data_dir);
    }
    if (tfield_enabled())
      io_tfd_.open("tfd" + sfx, prm.data_dir);
  }

  template <typename F>
  void operator()(const Grid_t& grid, F&& get_item)
  {
    static int pr_eval, pr_accum, pr_pfd, pr_tfd;
    if (!pr_eval) {
      pr_eval = prof_register("outf eval", 1., 0, 0);
      pr_accum = prof_register("outf accum", 1., 0, 0);
      pr_pfd = prof_register("outf pfd", 1., 0, 0);
      pr_tfd = prof_register("outf tfd", 1., 0, 0);
    }

    int timestep = grid.timestep();
    bool restarting_from_checkpoint = first_time_ && timestep != 0;
    if (restarting_from_checkpoint) {
      pfield_next_ = timestep + pfield_interval;

      if (tfield_enabled()) {
        // not counting initial tfd when timestep == tfield_first
        int n_tfds_already_written =
          (timestep - tfield_first) / tfield_interval;
        tfield_next_ =
          tfield_first + tfield_interval * (n_tfds_already_written + 1);
      }
    }
    first_time_ = false;

    bool do_pfield = pfield_interval > 0 && timestep >= pfield_next_;

    bool on_tfield_out_step = timestep >= tfield_next_;
    bool do_tfield = tfield_enabled() && on_tfield_out_step;

    bool in_tfield_averaging_range =
      tfield_next_ - timestep < tfield_average_length;
    bool on_tfield_averaging_step =
      (tfield_next_ - timestep) % tfield_average_every == 0;
    bool doaccum_tfield =
      tfield_enabled() && in_tfield_averaging_range && on_tfield_averaging_step;

    if (do_pfield || doaccum_tfield) {
      prof_start(pr_eval);
      auto&& item = get_item();
      auto&& pfd = psc::mflds::interior(grid, item.gt);
      prof_stop(pr_eval);

      if (do_pfield) {
        prof_start(pr_pfd);
        mpi_printf(grid.comm(), "***** Writing PFD output for '%s'\n",
                   item.name.c_str());
        pfield_next_ += pfield_interval;
        io_pfd_.write_step(grid, rn, rx, pfd, item.name, item.comp_names);
        prof_stop(pr_pfd);
      }

      if (doaccum_tfield) {
        if (!tfd_) {
          tfd_.reset(new Mfields{grid, item.gt.shape(3), {}});
        }

        prof_start(pr_accum);
        // tfd += pfd
        tfd_->gt() = tfd_->gt() + pfd;
        naccum_++;
        prof_stop(pr_accum);
      }

      if (do_tfield) {
        prof_start(pr_tfd);
        mpi_printf(grid.comm(), "***** Writing TFD output for '%s'\n",
                   item.name.c_str());
        tfield_next_ += tfield_interval;

        // convert accumulated values to correct temporal mean
        tfd_->gt() = (1. / naccum_) * tfd_->gt();

        // io_tfd_.begin_step(grid);
        // io_tfd_.set_subset(grid, rn, rx);
        // io_tfd_.write(tfd_->gt(), grid, item.name(),
        // item.comp_names()); io_tfd_.end_step();
        io_tfd_.write_step(grid, rn, rx, tfd_->gt(), item.name,
                           item.comp_names);
        naccum_ = 0;

        tfd_->gt().view() = 0;
        prof_stop(pr_tfd);
      }
    }
  }

private:
  int pfield_next_;
  int tfield_next_;
  Writer io_pfd_;
  Writer io_tfd_;
  std::unique_ptr<Mfields> tfd_;
  int naccum_ = 0;
  bool first_time_ =
    true; // to keep track so we can skip first output on restart
};

// ======================================================================
// OutputFieldsParams

struct OutputFieldsParams
{
  OutputFieldsItemParams fields;
  OutputFieldsItemParams moments;
};

// ======================================================================
// OutputFields

template <typename MfieldsState, typename Mparticles, typename Dim,
          typename Writer = WriterDefault>
class OutputFields
{
public:
  // ----------------------------------------------------------------------
  // ctor

  OutputFields(const Grid_t& grid, const OutputFieldsParams& prm)
    : fields{grid, prm.fields, ""}, moments{grid, prm.moments, "_moments"}
  {}

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MfieldsState& mflds, Mparticles& mprts)
  {
    const auto& grid = mflds._grid();

    static int pr, pr_fields, pr_moments;
    if (!pr) {
      pr = prof_register("outf", 1., 0, 0);
      pr_fields = prof_register("outf_fields", 1., 0, 0);
      pr_moments = prof_register("outf_moments", 1., 0, 0);
    }

    prof_start(pr);

    prof_start(pr_fields);
    fields(grid, [&]() {
      auto item = Item_jeh<MfieldsState>{};
      return make_mfields_gt(item(mflds), item.name(), item.comp_names());
    });
    prof_stop(pr_fields);

    prof_start(pr_moments);
    moments(grid, [&]() {
      auto item = Item_Moments<typename MfieldsState::Storage, Dim>(grid);
      return make_mfields_gt(item(mprts), item.name(), item.comp_names());
    });
    prof_stop(pr_moments);

    prof_stop(pr);
  };

public:
  OutputFieldsItem<Mfields_from_gt_t<Item_jeh<MfieldsState>>, Writer> fields;
  OutputFieldsItem<
    Mfields_from_gt_t<Item_Moments<typename MfieldsState::Storage, Dim>>,
    Writer>
    moments;
};
