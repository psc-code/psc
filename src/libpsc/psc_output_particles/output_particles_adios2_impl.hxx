#include "grid.hxx"
#include "output_particles.hxx"
#include "diagnostic_base.hxx"
#include "writer_adios2.hxx"

template <typename Mparticles>
class OutputParticlesAdios2
  : OutputParticlesBase
  , public ParticleDiagnosticBase<Mparticles>
{
  static OutputParticlesParams adjust_params(
    const OutputParticlesParams& params_in, const Grid_t& grid)
  {
    OutputParticlesParams params = params_in;
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
  OutputParticlesAdios2(const Grid_t& grid, const OutputParticlesParams& params)
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
  const OutputParticlesParams params_;
  WriterADIOS2 io_;
};
