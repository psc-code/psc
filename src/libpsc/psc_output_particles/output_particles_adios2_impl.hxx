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

template <typename Mparticles>
class OutputParticlesAdios2
  : OutputParticlesBase
  , public ParticleDiagnosticBase<Mparticles>
{
  using real_t = typename Mparticles::real_t;

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

    unsigned long local_n = 0; // gets set later

    if (params_.write_x)
      io_.DefineVariable<real_t>("x", {adios2::JoinedDim}, {}, {local_n});

    init_ = true;
  }

  void perform_diagnostic(Mparticles& mprts) override
  {
    const Grid_t& grid = mprts.grid();

    if (!init_) {
      init(grid);
    }

    //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\.
    // count number of particles of each kind
    int n_kinds = grid.kinds.size();
    std::vector<unsigned long> local_n_prts_per_kind(n_kinds);

    for (int p = 0; p < grid.n_patches(); p++) {
      for (auto prt = mprts.begin(p); prt != mprts.end(p); prt++) {
        local_n_prts_per_kind[prt->kind]++;
      }
    }
    ///////////////////////////////////////////////

    for (int kind_idx = 0; kind_idx < n_kinds; kind_idx++) {
      Grid_t::Kind kind = grid.kinds[kind_idx];
      unsigned long local_n_of_kind = local_n_prts_per_kind[kind_idx];

      std::string file_name =
        get_file_name(params_.basename, kind.name, grid.timestep());
      auto engine = io_.Open(file_name, adios2::Mode::Write);
      engine.BeginStep(adios2::StepMode::Append);

      //\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\.
      // data must live until EndStep(), so must be scoped accordingly

      std::vector<real_t> xs;

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // reserve as needed

      if (params_.write_x)
        xs.reserve(local_n_of_kind);

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // populate data

      for (int p = 0; p < grid.n_patches(); p++) {
        Grid_t::Patch patch = grid.patches[p];
        for (auto prt = mprts.begin(p); prt != mprts.end(p); prt++) {
          if (prt->kind != kind_idx)
            continue;

          if (params_.write_x)
            xs.push_back(prt->x[0] + patch.xb[0]);
        }
      }

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // schedule puts

      if (params_.write_x) {
        adios2::Variable<real_t> var = io_.InquireVariable<real_t>("x");
        var.SetSelection({{}, {local_n_of_kind}});
        engine.Put(var, xs.data(), adios2::Mode::Deferred);
      }

      //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      // performs puts - data must live until this point

      engine.EndStep();
      engine.Close();

      //////////////////////////////////////////////////////////////////////
    }
  }

  void operator()(Mparticles& mprts) { this->perform_diagnostic(mprts); }

private:
  const OutputParticlesAdios2Params params_;
  adios2::ADIOS adios_{MPI_COMM_WORLD};
  adios2::IO io_;
  bool init_ = false;
};
