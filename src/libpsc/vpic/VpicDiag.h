
#ifndef VPIC_DIAG_H
#define VPIC_DIAG_H

#include "vpic_iface.h"
#include "vpic.h"

#include "vpic/dumpmacros.h"

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

namespace dump_type {
  const int grid_dump = 0;
  const int field_dump = 1;
  const int hydro_dump = 2;
  const int particle_dump = 3;
  const int restart_dump = 4;
  const int history_dump = 5;
}

template<class FieldArray, class Particles, class Interpolator>
struct VpicDiagOps
{
  typedef VpicDiag Diag;

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
      
      /*--------------------------------------------------------------------------
       * Field data output
       *------------------------------------------------------------------------*/

      if(step == -1 || should_dump(fields)) field_dump(fa, simulation, diag.fdParams);
      
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
    
    if (!fname) ERROR(("Invalid file name"));
 
    if (rank==0) {
      status = fileIO.open(fname, append ? io_append : io_write);
      if (status==fail) ERROR(( "Could not open \"%s\".", fname));
      else {
	if (append==0) {
	  fileIO.print("%% Layout\n%% step ex ey ez bx by bz" );
	  LIST_FOR_EACH(sp, particles.sl_)
	    fileIO.print(" \"%s\"", sp->name);
	  fileIO.print("\n");
	  fileIO.print("%% timestep = %e\n", fa.g->dt);
	}
	fileIO.print( "%li", (long)simulation->step() );
      }
    }

    fa.energy_f(en_f);
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

  // ----------------------------------------------------------------------
  // field_dump
  
  void field_dump(FieldArray& fa, vpic_simulation* simulation, DumpParameters& dumpParams)
  {
#define step simulation->step
#define rank simulation->rank
#define nproc simulation->nproc
#define nxout simulation->nxout
#define nyout simulation->nyout
#define nzout simulation->nzout
#define dxout simulation->dxout
#define dyout simulation->dyout
#define dzout simulation->dzout
    grid_t* grid = simulation->grid;
    field_array_t* field_array = &fa;
    
    // Create directory for this time step
    char timeDir[256];
    sprintf(timeDir, "%s/T.%ld", dumpParams.baseDir, (long)step());
    simulation->dump_mkdir(timeDir);
  
    // Open the file for output
    char filename[256];
    sprintf(filename, "%s/T.%ld/%s.%ld.%d", dumpParams.baseDir, (long)step(),
	    dumpParams.baseFileName, (long)step(), rank());
  
    FileIO fileIO;
    FileIOStatus status;
  
    status = fileIO.open(filename, io_write);
    if( status==fail ) ERROR(( "Failed opening file: %s", filename ));
  
    // convenience
    const size_t istride(dumpParams.stride_x);
    const size_t jstride(dumpParams.stride_y);
    const size_t kstride(dumpParams.stride_z);
  
    // Check stride values.
    if(remainder(grid->nx, istride) != 0)
      ERROR(("x stride must be an integer factor of nx"));
    if(remainder(grid->ny, jstride) != 0)
      ERROR(("y stride must be an integer factor of ny"));
    if(remainder(grid->nz, kstride) != 0)
      ERROR(("z stride must be an integer factor of nz"));
  
    int dim[3];
  
    /* define to do C-style indexing */
# define f(x,y,z) f[ VOXEL(x,y,z, grid->nx,grid->ny,grid->nz) ]

    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    nxout = (grid->nx)/istride;
    nyout = (grid->ny)/jstride;
    nzout = (grid->nz)/kstride;
    dxout = (grid->dx)*istride;
    dyout = (grid->dy)*jstride;
    dzout = (grid->dz)*kstride;

    /* Banded output will write data as a single block-array as opposed to
     * the Array-of-Structure format that is used for native storage.
     *
     * Additionally, the user can specify a stride pattern to reduce
     * the resolution of the data that are output.  If a stride is
     * specified for a particular dimension, VPIC will write the boundary
     * plus every "stride" elements in that dimension. */

    if(dumpParams.format == band) {
    
      WRITE_HEADER_V0(dump_type::field_dump, -1, 0, fileIO);
    
      dim[0] = nxout+2;
      dim[1] = nyout+2;
      dim[2] = nzout+2;
    
      WRITE_ARRAY_HEADER(simulation->field_array->f, 3, dim, fileIO);
    
      // Create a variable list of field values to output.
      size_t numvars = std::min(dumpParams.output_vars.bitsum(),
				total_field_variables);
      size_t * varlist = new size_t[numvars];
    
      for(size_t i(0), c(0); i<total_field_variables; i++)
	if(dumpParams.output_vars.bitset(i)) varlist[c++] = i;
    
      // more efficient for standard case
      if(istride == 1 && jstride == 1 && kstride == 1)
	for(size_t v(0); v<numvars; v++) {
	  for(size_t k(0); k<nzout+2; k++) {
	    for(size_t j(0); j<nyout+2; j++) {
	      for(size_t i(0); i<nxout+2; i++) {
		const uint32_t * fref = reinterpret_cast<uint32_t *>(&field_array->f(i,j,k));
		fileIO.write(&fref[varlist[v]], 1);
	      }
	    }
	  }
	}

      else

	for(size_t v(0); v<numvars; v++) {
	  for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
	    for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
	      for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
		const uint32_t * fref = reinterpret_cast<uint32_t *>(&field_array->f(ioff,joff,koff));
		fileIO.write(&fref[varlist[v]], 1);
	      }
	    }
	  }
	}
    
      delete[] varlist;

    } else { // band_interleave
    
      WRITE_HEADER_V0(dump_type::field_dump, -1, 0, fileIO);
    
      dim[0] = nxout+2;
      dim[1] = nyout+2;
      dim[2] = nzout+2;
    
      WRITE_ARRAY_HEADER(field_array->f, 3, dim, fileIO);
    
      if(istride == 1 && jstride == 1 && kstride == 1)
	fileIO.write(field_array->f, dim[0]*dim[1]*dim[2]);
      else
	for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
	  for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
	    for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
	      fileIO.write(&field_array->f(ioff,joff,koff), 1);
	    }
	  }
	}
    }
  
# undef f
    
    if( fileIO.close() ) ERROR(( "File close failed on field dump!!!" ));
#undef step
#undef rank
#undef nproc
#undef nxout
#undef nyout
#undef nzout
  }
  
};


#endif
