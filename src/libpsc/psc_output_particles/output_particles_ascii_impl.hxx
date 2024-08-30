
#include "output_particles.hxx"

struct OutputParticlesAscii
  : OutputParticlesParams
  , OutputParticlesBase
{
  OutputParticlesAscii(const Grid_t& grid, const OutputParticlesParams& params)
    : OutputParticlesParams(params), comm_{grid.comm()}
  {}

  // ----------------------------------------------------------------------
  // operator()

  template <typename Mparticles>
  void operator()(Mparticles& mprts)
  {
    const auto& grid = mprts.grid();
    if (every_step < 0 || grid.timestep() % every_step != 0) {
      return;
    }

    int rank;
    MPI_Comm_rank(comm_, &rank);
    int slen = strlen(data_dir) + strlen(basename) + 21;
    char filename[slen];
    snprintf(filename, slen, "%s/%s.%06d_p%06d.asc", data_dir, basename,
             grid.timestep(), rank);

    FILE* file = fopen(filename, "w");
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      int n = 0;
      for (auto prt : accessor[p]) {
        fprintf(file, "%d %g %g %g %g %g %g %g %d\n", n, prt.x()[0], prt.x()[1],
                prt.x()[2], prt.u()[0], prt.u()[1], prt.u()[2], prt.w(),
                prt.kind());
        n++;
      }
    }

    fclose(file);
  }

private:
  MPI_Comm comm_;
};
