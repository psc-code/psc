
#include "psc_particles_double.h"
#include "output_particles.hxx"

#include <mrc_params.h>
#include <string.h>

#include <string.h>

#define to_psc_output_particles_ascii(out) \
  mrc_to_subobj(out, struct psc_output_particles_ascii)

struct psc_output_particles_ascii : OutputParticlesParams, OutputParticlesBase
{
  psc_output_particles_ascii(const Grid_t& grid, const OutputParticlesParams& params)
    : OutputParticlesParams(params),
      comm_{grid.comm()}
  {}

  // ----------------------------------------------------------------------
  // run

  void run(MparticlesBase& mprts_base) override
  {
    const auto& grid = mprts_base.grid();
    if (every_step < 0 || grid.timestep() % every_step != 0) {
      return;
    }
    
    int rank;
    MPI_Comm_rank(comm_, &rank);
    char filename[strlen(data_dir) + strlen(basename) + 19];
    sprintf(filename, "%s/%s.%06d_p%06d.asc", data_dir,
	    basename, grid.timestep(), rank);
    
    auto& mprts = mprts_base.get_as<MparticlesDouble>();
    
    FILE *file = fopen(filename, "w");
    for (int p = 0; p < mprts.n_patches(); p++) {
      int n = 0;
      for (auto& prt : mprts[p]) {
	fprintf(file, "%d %g %g %g %g %g %g %g %d\n",
		n, prt.x[0], prt.x[1], prt.x[2],
		prt.p[0], prt.p[1], prt.p[2],
		prt.w, prt.kind);
	n++;
      }
    }
    
    mprts_base.put_as(mprts, MP_DONT_COPY);
    
    fclose(file);
  }

private:
  MPI_Comm comm_;
};

// ======================================================================
// psc_output_particles: subclass "ascii"

struct psc_output_particles_ops_ascii : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_ascii>;
  psc_output_particles_ops_ascii() {
    name                  = "ascii";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
  }
} psc_output_particles_ascii_ops;
