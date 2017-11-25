
#ifndef VPIC_DIAG_H
#define VPIC_DIAG_H

#include "vpic_iface.h"
#include "vpic.h"

// ----------------------------------------------------------------------
// VpicDiag

struct VpicDiag {
  int interval;
  int energies_interval;
  int fields_interval;
  int ehydro_interval;
  int Hhydro_interval;
  int eparticle_interval;
  int Hparticle_interval;
  int restart_interval;

  // state
  int rtoggle;               // enables save of last 2 restart dumps for safety
  // Output variables
  DumpParameters fdParams;
  DumpParameters hedParams;
  DumpParameters hHdParams;
  std::vector<DumpParameters *> outputParams;

  VpicDiag(vpic_simulation* simulation);
  void init(int interval_);
  void setup();
  void run();

private:
  vpic_simulation* simulation_;
};



void vpic_simulation_diagnostics(vpic_simulation *simulation, VpicDiag *diag);
void vpic_simulation_setup_diagnostics(vpic_simulation *simulation, VpicDiag *diag);


inline VpicDiag::VpicDiag(vpic_simulation *simulation)
  : simulation_(simulation)
{
}

inline void VpicDiag::init(int interval_)
{
  rtoggle = 0;
  
  interval = interval_;
  fields_interval = interval_;
  ehydro_interval = interval_;
  Hhydro_interval = interval_;
  eparticle_interval = 8 * interval_;
  Hparticle_interval = 8 * interval_;
  restart_interval = 8000;

  energies_interval = 50;

  MPI_Comm comm = MPI_COMM_WORLD;
  mpi_printf(comm, "interval = %d\n", interval);
  mpi_printf(comm, "energies_interval: %d\n", energies_interval);
}

inline void VpicDiag::setup()
{
  vpic_simulation_setup_diagnostics(simulation_, this);
}

inline void VpicDiag::run()
{
}

// ======================================================================
// VpicDiagOps

template<class FieldArrayOps, class ParticlesOps, class InterpolatorOps>
struct VpicDiagOps
{
  typedef VpicDiag Diag;
  typedef typename FieldArrayOps::FieldArray FieldArray;
  typedef typename ParticlesOps::Particles Particles;
  typedef typename InterpolatorOps::Interpolator Interpolator;

#define should_dump(x)                                                  \
  (diag.x##_interval>0 && remainder(step, diag.x##_interval) == 0)

  void diagnostics_run(VpicDiag& diag, FieldArray& fa, Particles& particles,
		       Interpolator& interpolator)
  {
    TIC {
      vpic_simulation *simulation = diag.simulation_;
      int64_t step = simulation->step();
      
      /*--------------------------------------------------------------------------
       * Data output directories
       * WARNING: The directory list passed to "global_header" must be
       * consistent with the actual directories where fields and species are
       * output using "field_dump" and "hydro_dump".
       *
       * DIRECTORY PATHES SHOULD BE RELATIVE TO
       * THE LOCATION OF THE GLOBAL HEADER!!!
       *------------------------------------------------------------------------*/
      
      // Normal rundata dump
      if (step==0) {
	simulation->dump_mkdir("fields");
	simulation->dump_mkdir("hydro");
	simulation->dump_mkdir("rundata");
	simulation->dump_mkdir("injectors");
	simulation->dump_mkdir("restart1");  // 1st backup
	simulation->dump_mkdir("restart2");  // 2nd backup
	simulation->dump_mkdir("particle");
	
	simulation->dump_grid("rundata/grid");
	simulation->dump_materials("rundata/materials");
	simulation->dump_species("rundata/species");
	simulation->global_header("global", diag.outputParams);
      }

      // Normal rundata energies dump
      if(should_dump(energies)) {
	dump_energies(simulation, "rundata/energies", step == 0 ? 0 : 1,
		      fa, particles, interpolator);
      }
      
      vpic_simulation_diagnostics(simulation, &diag);
    } TOC(user_diagnostics, 1);
  }

#undef should_dump

  // ----------------------------------------------------------------------
  // dump_energies

  void dump_energies(vpic_simulation *simulation, const char *fname, int append,
		     FieldArray& fa, Particles& particles,
		     Interpolator& interpolator)
  {
    double en_f[6], en_p;
    species_t *sp;
    FileIO fileIO;
    FileIOStatus status(fail);

    int rank = simulation->rank();
    
    if( !fname ) ERROR(("Invalid file name"));
 
    if( rank==0 ) {
      status = fileIO.open(fname, append ? io_append : io_write);
      if(status==fail) ERROR(( "Could not open \"%s\".", fname));
      else {
	if( append==0 ) {
	  fileIO.print("%% Layout\n%% step ex ey ez bx by bz" );
	  LIST_FOR_EACH(sp, particles.sl_)
	    fileIO.print(" \"%s\"", sp->name);
	  fileIO.print("\n");
	  fileIO.print("%% timestep = %e\n", fa.g->dt);
	}
	fileIO.print( "%li", (long)simulation->step() );
      }
    }
 
    fa.kernel->energy_f(en_f, &fa);
    if (rank==0 && status!=fail )
      fileIO.print( " %e %e %e %e %e %e",
		    en_f[0], en_f[1], en_f[2],
		    en_f[3], en_f[4], en_f[5] );
 
    LIST_FOR_EACH(sp, particles.sl_) {
      en_p = energy_p(sp, &interpolator);
      if (rank==0 && status != fail) fileIO.print(" %e", en_p);
    }
 
    if (rank==0 && status != fail) {
      fileIO.print("\n");
      if (fileIO.close()) ERROR(("File close failed on dump energies!!!"));
    }
  }
 
};

#endif
