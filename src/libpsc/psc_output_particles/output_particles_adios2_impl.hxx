#include "grid.hxx"
#include "output_particles.hxx"
#include "diagnostic_base.hxx"
#include "writer_adios2.hxx"

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

  void perform_diagnostic(Mparticles& mprts) override
  {
    if (!io_) {
      io_.open(params_.basename, params_.data_dir);
    }
    // TODO
  }

  void operator()(Mparticles& mprts) { this->perform_diagnostic(mprts); }

private:
  const OutputParticlesAdios2Params params_;
  WriterADIOS2 io_;
};
