
#include "output_particles.hxx"

#include "psc_particles_double.h"

struct OutputParticlesAscii : OutputParticlesParams, OutputParticlesBase
{
  OutputParticlesAscii(const Grid_t& grid, const OutputParticlesParams& params)
    : OutputParticlesParams(params),
      comm_{grid.comm()}
  {}

  // ----------------------------------------------------------------------
  // operator()

  void operator()(MparticlesBase& mprts_base) 
  {
    const auto& grid = mprts_base.grid();
    if (every_step < 0 || grid.timestep() % every_step != 0) {
      return;
    }
    
    int rank;
    MPI_Comm_rank(comm_, &rank);
    char filename[strlen(data_dir) + strlen(basename) + 21];
    sprintf(filename, "%s/%s.%06d_p%06d.asc", data_dir,
	    basename, grid.timestep(), rank);
    
    auto& mprts = mprts_base.get_as<MparticlesDouble>();
    
    FILE *file = fopen(filename, "w");
    auto accessor = mprts.accessor();
    for (int p = 0; p < mprts.n_patches(); p++) {
      int n = 0;
      for (auto prt : accessor[p]) {
	fprintf(file, "%d %g %g %g %g %g %g %g %d\n",
		n, prt.x()[0], prt.x()[1], prt.x()[2],
		prt.u()[0], prt.u()[1], prt.u()[2],
		prt.w(), prt.kind()); 
	n++;
      }
    }
    
    mprts_base.put_as(mprts, MP_DONT_COPY);
    
    fclose(file);
  }

private:
  MPI_Comm comm_;
};

