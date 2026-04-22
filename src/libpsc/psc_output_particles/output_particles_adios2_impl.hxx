#include <adios2.h>

#include "grid.hxx"
#include "output_particles.hxx"
#include "diagnostic_base.hxx"

struct OutputParticlesAdios2Params : OutputParticlesParams
{
  bool write_x = true;
  bool write_y = true;
  bool write_z = true;
  bool write_px = true;
  bool write_py = true;
  bool write_pz = true;
  bool write_q = false;
  bool write_m = false;
  bool write_w = true;
  bool write_id = false;
  bool write_tag = false;

  OutputParticlesAdios2Params() {}

  OutputParticlesAdios2Params(OutputParticlesParams params)
    : OutputParticlesParams{params}
  {}
};

std::string get_file_name(std::string basename, std::string kind_name, int step)
{
  int min_step_digits = 9;
  std::string padded_step = std::to_string(step);
  if (padded_step.length() < min_step_digits) {
    padded_step.insert(0, min_step_digits - padded_step.length(), '0');
  }
  return basename + "." + kind_name + "." + padded_step + ".bp";
}

template <typename Mparticles,
          typename RealOverride = typename Mparticles::real_t>
class OutputParticlesAdios2
  : OutputParticlesBase
  , public ParticleDiagnosticBase<Mparticles>
{
  using real_t = RealOverride;

  static OutputParticlesAdios2Params adjust_params(
    const OutputParticlesAdios2Params& params_in, const Grid_t& grid)
  {
    OutputParticlesAdios2Params params = params_in;
    for (int d = 0; d < 3; d++) {
      if (params.hi[d] == 0) {
        params.hi[d] = grid.domain.gdims[d];
      }
      assert(params.lo[d] >= 0);
      assert(params.hi[d] <= grid.domain.gdims[d]);
    }
    return params;
  }

public:
  OutputParticlesAdios2(const Grid_t& grid,
                        const OutputParticlesAdios2Params& params)
    : params_{adjust_params(params, grid)}
  {}

  void init(const Grid_t& grid)
  {
    io_ = adios_.DeclareIO("PrtWriter");
    // TODO set IO parameters, configured by params_
    // (the default options are generally better at large scales, but doesn't
    // take advantage of relatively higher file bandwidth at smaller scales)

    unsigned long local_n = 0; // gets set later

    if (params_.write_x)
      io_.DefineVariable<real_t>("x", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_y)
      io_.DefineVariable<real_t>("y", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_z)
      io_.DefineVariable<real_t>("z", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_px)
      io_.DefineVariable<real_t>("px", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_py)
      io_.DefineVariable<real_t>("py", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_pz)
      io_.DefineVariable<real_t>("pz", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_q)
      io_.DefineVariable<float>("q", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_m)
      io_.DefineVariable<float>("m", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_w)
      io_.DefineVariable<float>("w", {adios2::JoinedDim}, {}, {local_n});
    if (params_.write_id)
      io_.DefineVariable<psc::particle::Id>("id", {adios2::JoinedDim}, {},
                                            {local_n});
    if (params_.write_tag)
      io_.DefineVariable<psc::particle::Tag>("tag", {adios2::JoinedDim}, {},
                                             {local_n});

    io_.DefineAttribute("gdims", grid.domain.gdims.data(),
                        grid.domain.gdims.size());
    io_.DefineAttribute("corner", grid.domain.corner.data(),
                        grid.domain.corner.size());
    io_.DefineAttribute("length", grid.domain.length.data(),
                        grid.domain.length.size());

    bool allow_modification{true};
    io_.DefineAttribute("step", grid.timestep(), "", "/", allow_modification);
    io_.DefineAttribute("time", grid.timestep() * grid.dt, "", "/",
                        allow_modification);

    init_ = true;
  }

  void perform_diagnostic(Mparticles& mprts) override
  {
    static int pr_all, pr_count, pr_fill, pr_schedule, pr_write;
    if (!pr_all) {
      pr_all = prof_register("outp_adios2", 1., 0, 0);
      pr_count = prof_register("outp_adios2: count", 1., 0, 0);
      pr_fill = prof_register("outp_adios2: fill", 1., 0, 0);
      pr_schedule = prof_register("outp_adios2: schedule", 1., 0, 0);
      pr_write = prof_register("outp_adios2: write", 1., 0, 0);
    }

    const Grid_t& grid = mprts.grid();

    if (params_.every_step <= 0 || grid.timestep() % params_.every_step != 0) {
      return;
    }

    prof_start(pr_all);

    if (!init_) {
      init(grid);
    }

    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\.
    // count number of particles of each kind
    prof_start(pr_count);

    int n_kinds = grid.kinds.size();
    std::vector<unsigned long> local_n_prts_per_kind(n_kinds);

    for (int p = 0; p < grid.n_patches(); p++) {
      for (auto prt = mprts.begin(p); prt != mprts.end(p); prt++) {
        local_n_prts_per_kind[prt->kind]++;
      }
    }

    prof_stop(pr_count);
    ///////////////////////////////////////////////

    for (int kind_idx = 0; kind_idx < n_kinds; kind_idx++) {
      Grid_t::Kind kind = grid.kinds[kind_idx];
      unsigned long local_n_of_kind = local_n_prts_per_kind[kind_idx];

      mpi_printf(grid.comm(), "***** Writing PRT output for '%s'\n",
                 kind.name.c_str());

      std::string file_name =
        get_file_name(params_.basename, kind.name, grid.timestep());
      auto engine = io_.Open(file_name, adios2::Mode::Write);
      engine.BeginStep(adios2::StepMode::Append);

      //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\.
      // data must live until EndStep(), so must be scoped accordingly

      std::vector<real_t> xs;
      std::vector<real_t> ys;
      std::vector<real_t> zs;
      std::vector<real_t> pxs;
      std::vector<real_t> pys;
      std::vector<real_t> pzs;
      std::vector<float> qs;
      std::vector<float> ms;
      std::vector<float> ws;
      std::vector<psc::particle::Id> ids;
      std::vector<psc::particle::Tag> tags;

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // reserve as needed
      if (kind_idx == 0)
        prof_start(pr_fill);
      else
        prof_restart(pr_fill);

      // TODO OPT maybe just keep these temps around instead of reallocating
      // on every write

      if (params_.write_x)
        xs.reserve(local_n_of_kind);
      if (params_.write_y)
        ys.reserve(local_n_of_kind);
      if (params_.write_z)
        zs.reserve(local_n_of_kind);
      if (params_.write_px)
        pxs.reserve(local_n_of_kind);
      if (params_.write_py)
        pys.reserve(local_n_of_kind);
      if (params_.write_pz)
        pzs.reserve(local_n_of_kind);
      if (params_.write_q)
        qs.reserve(local_n_of_kind);
      if (params_.write_m)
        ms.reserve(local_n_of_kind);
      if (params_.write_w)
        ws.reserve(local_n_of_kind);
      if (params_.write_id)
        ids.reserve(local_n_of_kind);
      if (params_.write_tag)
        tags.reserve(local_n_of_kind);

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // populate data

      for (int p = 0; p < grid.n_patches(); p++) {
        Grid_t::Patch patch = grid.patches[p];
        for (auto prt = mprts.begin(p); prt != mprts.end(p); prt++) {
          if (prt->kind != kind_idx)
            continue;

          // FIXME OPT: this is 11 ifs in a hot loop. Kinda just hoping that
          // branch prediction and/or cache hits make this ok.

          if (params_.write_x)
            xs.push_back(prt->x[0] + patch.xb[0]);
          if (params_.write_y)
            ys.push_back(prt->x[1] + patch.xb[1]);
          if (params_.write_z)
            zs.push_back(prt->x[2] + patch.xb[2]);
          if (params_.write_px)
            pxs.push_back(prt->u[0]);
          if (params_.write_py)
            pys.push_back(prt->u[1]);
          if (params_.write_pz)
            pzs.push_back(prt->u[2]);
          if (params_.write_q)
            qs.push_back(kind.q);
          if (params_.write_m)
            ms.push_back(kind.m);
          if (params_.write_w)
            ws.push_back(prt->qni_wni / kind.q);
          if (params_.write_id)
            ids.push_back(prt->id());
          if (params_.write_tag)
            tags.push_back(prt->tag());
        }
      }

      prof_stop(pr_fill);
      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // schedule puts
      if (kind_idx == 0)
        prof_start(pr_schedule);
      else
        prof_restart(pr_schedule);

      if (params_.write_x) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("x");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, xs.data(), adios2::Mode::Deferred);
      }
      if (params_.write_y) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("y");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, ys.data(), adios2::Mode::Deferred);
      }
      if (params_.write_z) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("z");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, zs.data(), adios2::Mode::Deferred);
      }
      if (params_.write_px) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("px");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, pxs.data(), adios2::Mode::Deferred);
      }
      if (params_.write_py) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("py");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, pys.data(), adios2::Mode::Deferred);
      }
      if (params_.write_pz) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("pz");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, pzs.data(), adios2::Mode::Deferred);
      }
      if (params_.write_q) {
        adios2::Variable<float> var = io_.InquireVariable<float>("q");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, qs.data(), adios2::Mode::Deferred);
      }
      if (params_.write_m) {
        adios2::Variable<float> var = io_.InquireVariable<float>("m");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, ms.data(), adios2::Mode::Deferred);
      }
      if (params_.write_w) {
        adios2::Variable<float> var = io_.InquireVariable<float>("w");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, ws.data(), adios2::Mode::Deferred);
      }
      if (params_.write_id) {
        adios2::Variable<psc::particle::Id> var =
          io_.InquireVariable<psc::particle::Id>("id");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, ids.data(), adios2::Mode::Deferred);
      }
      if (params_.write_tag) {
        adios2::Variable<psc::particle::Tag> var =
          io_.InquireVariable<psc::particle::Tag>("tag");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, tags.data(), adios2::Mode::Deferred);
      }

      prof_stop(pr_schedule);
      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // update io-level attributes

      io_.DefineAttribute("step", grid.timestep());
      io_.DefineAttribute("time", grid.timestep() * grid.dt);

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // performs puts - data must live until this point
      if (kind_idx == 0)
        prof_start(pr_write);
      else
        prof_restart(pr_write);

      engine.EndStep();
      engine.Close();

      prof_stop(pr_write);
      //////////////////////////////////////////////////////////////////////
    }

    prof_stop(pr_all);
  }

  void operator()(Mparticles& mprts) { this->perform_diagnostic(mprts); }

private:
  const OutputParticlesAdios2Params params_;
  adios2::ADIOS adios_{MPI_COMM_WORLD};
  adios2::IO io_;
  bool init_ = false;
};
