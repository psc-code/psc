
#pragma once

#include "diagnostic_base.hxx"
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
// BaseOutputFieldItemParams

struct BaseOutputFieldItemParams
{
  int out_interval = 0; // difference between output timesteps (0 = disable)
  std::string data_dir = ".";
  Int3 rn = {};
  Int3 rx = {10000000, 10000000, 10000000};

  bool enabled() { return out_interval > 0; }

  bool do_out(int timestep)
  {
    return enabled() && timestep % out_interval == 0;
  }
};

struct OutputPfieldItemParams : BaseOutputFieldItemParams
{};

struct OutputTfieldItemParams : BaseOutputFieldItemParams
{
  // max range of timesteps over which to average (capped at `out_interval`)
  int average_length = 1000000;
  // difference between timesteps used for average
  int sample_interval = 1;

  // Returns whether to accumulate on this timestep.
  bool do_accum(int timestep)
  {
    // next_out could be this timestep
    int n_intervals_elapsed = (timestep - 1) / out_interval;
    int next_out = out_interval * (n_intervals_elapsed + 1);

    bool in_averaging_range = next_out - timestep < average_length;
    bool on_averaging_step = (next_out - timestep) % sample_interval == 0;
    return enabled() && in_averaging_range && on_averaging_step;
  }
};

// ======================================================================
// GetItem*

struct GetItemJeh
{
  static std::string suffix() { return ""; }

  template <typename Mparticles, typename MfieldsState>
  static auto get_item(Mparticles& mprts, MfieldsState& mflds)
  {
    auto item = Item_jeh<MfieldsState>{};
    return make_mfields_gt(item(mflds), item.name(), item.comp_names());
  }
};

template <typename Dim>
struct GetItemMoments
{
  static std::string suffix() { return "_moments"; }

  template <typename Mparticles, typename MfieldsState>
  static auto get_item(Mparticles& mprts, MfieldsState& mflds)
  {
    auto item = Item_Moments<typename MfieldsState::Storage, Dim>(mprts.grid());
    return make_mfields_gt(item(mprts), item.name(), item.comp_names());
  }
};

// ======================================================================
// OutputFieldsItemParams

struct OutputFieldsItemParams
{
  OutputPfieldItemParams pfield;
  OutputTfieldItemParams tfield;
};

// ======================================================================
// OutputFieldsItem

// TODO infer Mfields from MfieldsState and/or GetItem
template <typename Mfields, typename MfieldsState, typename Mparticles,
          typename GetItem, typename Writer = WriterDefault>
class OutputFieldsItem
  : public OutputFieldsItemParams
  , public DiagnosticBase<Mparticles, MfieldsState>
{
public:
  OutputFieldsItem() = default;

  OutputFieldsItem(const OutputFieldsItemParams& prm)
    : OutputFieldsItemParams{prm}
  {}

  void perform_diagnostic(Mparticles& mprts, MfieldsState& mflds) override
  {
    const Grid_t& grid = mflds.grid();

    static int pr_outf, pr_eval, pr_accum, pr_pfd, pr_tfd;
    if (!pr_outf) {
      pr_outf = prof_register(("outf" + GetItem::suffix()).c_str(), 1., 0, 0);
      pr_eval = prof_register("outf eval", 1., 0, 0);
      pr_accum = prof_register("outf accum", 1., 0, 0);
      pr_pfd = prof_register("outf pfd", 1., 0, 0);
      pr_tfd = prof_register("outf tfd", 1., 0, 0);
    }

    prof_start(pr_outf);

    int timestep = grid.timestep();

    bool do_pfield = pfield.do_out(timestep);
    bool do_tfield = tfield.do_out(timestep);
    bool do_tfield_accum = tfield.do_accum(timestep);

    if (do_pfield || do_tfield_accum) {
      prof_start(pr_eval);
      auto&& item = GetItem::get_item(mprts, mflds);
      auto&& pfd = psc::mflds::interior(grid, item.gt);
      prof_stop(pr_eval);

      if (do_pfield) {
        if (!io_pfd_) {
          io_pfd_.open("pfd" + GetItem::suffix(), pfield.data_dir);
        }

        prof_start(pr_pfd);
        mpi_printf(grid.comm(), "***** Writing PFD output for '%s'\n",
                   item.name.c_str());
        io_pfd_.write_step(grid, pfield.rn, pfield.rx, pfd, item.name,
                           item.comp_names);
        prof_stop(pr_pfd);
      }

      if (do_tfield_accum) {
        if (!tfd_) {
          tfd_.reset(new Mfields{grid, item.gt.shape(3), {}});
        }

        prof_start(pr_accum);
        tfd_->gt() = tfd_->gt() + pfd;
        naccum_++;
        prof_stop(pr_accum);
      }

      if (do_tfield) {
        if (!io_tfd_) {
          io_tfd_.open("tfd" + GetItem::suffix(), tfield.data_dir);
        }

        prof_start(pr_tfd);
        mpi_printf(grid.comm(), "***** Writing TFD output for '%s'\n",
                   item.name.c_str());

        // convert accumulated values to correct temporal mean
        tfd_->gt() = (1. / naccum_) * tfd_->gt();

        io_tfd_.write_step(grid, tfield.rn, tfield.rx, tfd_->gt(), item.name,
                           item.comp_names);
        naccum_ = 0;

        tfd_->gt().view() = 0;
        prof_stop(pr_tfd);
      }
    }

    prof_stop(pr_outf);
  }

private:
  Writer io_pfd_;
  Writer io_tfd_;
  std::unique_ptr<Mfields> tfd_;
  int naccum_ = 0;
};

template <typename MfieldsState, typename Mparticles,
          typename Writer = WriterDefault>
using OutputFields =
  OutputFieldsItem<Mfields_from_gt_t<Item_jeh<MfieldsState>>, MfieldsState,
                   Mparticles, GetItemJeh, Writer>;

template <typename MfieldsState, typename Mparticles, typename Dim,
          typename Writer = WriterDefault>
using OutputMoments = OutputFieldsItem<
  Mfields_from_gt_t<Item_Moments<typename MfieldsState::Storage, Dim>>,
  MfieldsState, Mparticles, GetItemMoments<Dim>, Writer>;