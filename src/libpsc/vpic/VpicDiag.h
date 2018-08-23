
#ifndef VPIC_DIAG_H
#define VPIC_DIAG_H

#include "vpic/vpic.h"
#include "util/io/FileUtils.h"
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
};


// ======================================================================
// VpicDiagOps

#define _WRITE_HEADER_V0(dump_type,sp_id,q_m,fileIO) do { \
    /* Binary compatibility information */               \
    WRITE( char,      CHAR_BIT,               fileIO );  \
    WRITE( char,      sizeof(short int),      fileIO );  \
    WRITE( char,      sizeof(int),            fileIO );  \
    WRITE( char,      sizeof(float),          fileIO );  \
    WRITE( char,      sizeof(double),         fileIO );  \
    WRITE( short int, 0xcafe,                 fileIO );  \
    WRITE( int,       0xdeadbeef,             fileIO );  \
    WRITE( float,     1.0,                    fileIO );  \
    WRITE( double,    1.0,                    fileIO );  \
    /* Dump type and header format version */            \
    WRITE( int,       0 /* Version */,        fileIO );  \
    WRITE( int,       dump_type,              fileIO );  \
    /* High level information */                         \
    WRITE( int,       step,                   fileIO );  \
    WRITE( int,       nxout,                  fileIO );  \
    WRITE( int,       nyout,                  fileIO );  \
    WRITE( int,       nzout,                  fileIO );  \
    WRITE( float,     grid->dt,               fileIO );  \
    WRITE( float,     dxout,                  fileIO );  \
    WRITE( float,     dyout,                  fileIO );  \
    WRITE( float,     dzout,                  fileIO );  \
    WRITE( float,     grid->x0,               fileIO );  \
    WRITE( float,     grid->y0,               fileIO );  \
    WRITE( float,     grid->z0,               fileIO );  \
    WRITE( float,     grid->cvac,             fileIO );  \
    WRITE( float,     grid->eps0,             fileIO );  \
    WRITE( float,     0 /* damp */,           fileIO );  \
    WRITE( int,       rank,                   fileIO );  \
    WRITE( int,       nproc,                  fileIO );  \
    /* Species parameters */                             \
    WRITE( int,       sp_id,                  fileIO );  \
    WRITE( float,     q_m,                    fileIO );  \
  } while(0)
 
namespace dump_type {
  const int grid_dump = 0;
  const int field_dump = 1;
  const int hydro_dump = 2;
  const int particle_dump = 3;
  const int restart_dump = 4;
  const int history_dump = 5;
}

static FieldInfo fieldInfo[12] = {
	{ "Electric Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Electric Field Divergence Error", "SCALAR", "1", "FLOATING_POINT",
		sizeof(float) },
	{ "Magnetic Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Magnetic Field Divergence Error", "SCALAR", "1", "FLOATING_POINT",
		sizeof(float) },
	{ "TCA Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Bound Charge Density", "SCALAR", "1", "FLOATING_POINT", sizeof(float) },
	{ "Free Current Field", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Charge Density", "SCALAR", "1", "FLOATING_POINT", sizeof(float) },
	{ "Edge Material", "VECTOR", "3", "INTEGER", sizeof(material_id) },
	{ "Node Material", "SCALAR", "1", "INTEGER", sizeof(material_id) },
	{ "Face Material", "VECTOR", "3", "INTEGER", sizeof(material_id) },
	{ "Cell Material", "SCALAR", "1", "INTEGER", sizeof(material_id) }
}; // fieldInfo

static HydroInfo hydroInfo[5] = {
	{ "Current Density", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Charge Density", "SCALAR", "1", "FLOATING_POINT", sizeof(float) },
	{ "Momentum Density", "VECTOR", "3", "FLOATING_POINT", sizeof(float) },
	{ "Kinetic Energy Density", "SCALAR", "1", "FLOATING_POINT",
		sizeof(float) },
	{ "Stress Tensor", "TENSOR", "6", "FLOATING_POINT", sizeof(float) }
}; // hydroInfo

#undef sim_log
#define sim_log(x) do {                                \
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);    \
    if( rank==0 ) {				       \
      std::cerr << "SIM_LOG: " << x << std::endl;      \
      std::cerr.flush();                               \
    }                                                  \
  } while(0)

template<class Mparticles, class MfieldsState, class MfieldsInterpolator, class MfieldsHydro,
	 class DiagOps, class ParticlesOps, class HydroArrayOps>
struct VpicDiagMixin
{
  using Particles = typename Mparticles::Particles;
  using Grid = typename Particles::Grid;
  
  void diagnostics_init(int interval_)
  {
    diag_.rtoggle = 0;
    
    diag_.interval = interval_;
    diag_.fields_interval = interval_;
    diag_.ehydro_interval = interval_;
    diag_.Hhydro_interval = interval_;
    diag_.eparticle_interval = 8 * interval_;
    diag_.Hparticle_interval = 8 * interval_;
    diag_.restart_interval = 8000;
    
    diag_.energies_interval = 50;
    
    MPI_Comm comm = MPI_COMM_WORLD;
    mpi_printf(comm, "interval = %d\n", diag_.interval);
    mpi_printf(comm, "energies_interval: %d\n", diag_.energies_interval);
  }

  void diagnostics_setup()
  {
    diag_.fdParams.format = band;
    sim_log( "Fields output format = band" );
    
    diag_.hedParams.format = band;
    sim_log( "Electron species output format = band" );
    
    diag_.hHdParams.format = band;
    sim_log( "Ion species output format = band" );
    
    // relative path to fields data from global header
    sprintf(diag_.fdParams.baseDir, "fields");

    // base file name for fields output
    sprintf(diag_.fdParams.baseFileName, "fields");

    diag_.fdParams.stride_x = 1;
    diag_.fdParams.stride_y = 1;
    diag_.fdParams.stride_z = 1;

    // add field parameters to list
    diag_.outputParams.push_back(&diag_.fdParams);

    sim_log( "Fields x-stride " << diag_.fdParams.stride_x );
    sim_log( "Fields y-stride " << diag_.fdParams.stride_y );
    sim_log( "Fields z-stride " << diag_.fdParams.stride_z );

    // relative path to electron species data from global header
    sprintf(diag_.hedParams.baseDir, "hydro");

    // base file name for fields output
    sprintf(diag_.hedParams.baseFileName, "ehydro");

    diag_.hedParams.stride_x = 1;
    diag_.hedParams.stride_y = 1;
    diag_.hedParams.stride_z = 1;

    // add electron species parameters to list
    diag_.outputParams.push_back(&diag_.hedParams);

    sim_log( "Electron species x-stride " << diag_.hedParams.stride_x );
    sim_log( "Electron species y-stride " << diag_.hedParams.stride_y );
    sim_log( "Electron species z-stride " << diag_.hedParams.stride_z );

    // relative path to ion species data from global header
    sprintf(diag_.hHdParams.baseDir, "hydro");

    // base file name for fields output
    sprintf(diag_.hHdParams.baseFileName, "Hhydro");

    diag_.hHdParams.stride_x = 1;
    diag_.hHdParams.stride_y = 1;
    diag_.hHdParams.stride_z = 1;

    sim_log( "Ion species x-stride " << diag_.hHdParams.stride_x );
    sim_log( "Ion species y-stride " << diag_.hHdParams.stride_y );
    sim_log( "Ion species z-stride " << diag_.hHdParams.stride_z );

    // add ion species parameters to list
    diag_.outputParams.push_back(&diag_.hHdParams);

    diag_.fdParams.output_variables( all );
    diag_.hedParams.output_variables( all );
    diag_.hHdParams.output_variables( all );

    char varlist[512];
    create_field_list(varlist, diag_.fdParams);

    sim_log( "Fields variable list: " << varlist );

    create_hydro_list(varlist, diag_.hedParams);

    sim_log( "Electron species variable list: " << varlist );

    create_hydro_list(varlist, diag_.hHdParams);

    sim_log( "Ion species variable list: " << varlist );
  }

#define should_dump(x)                                                  \
  (diag_.x##_interval>0 && remainder(step, diag_.x##_interval) == 0)

  void diagnostics_run(Mparticles& mprts, MfieldsState& mflds,
		       MfieldsInterpolator& interpolator, MfieldsHydro& mflds_hydro, const int np[3])
  {
    TIC {
      const Grid* g = mflds.vgrid();
      int64_t step = g->step;
      
      // Normal rundata dump
      if (step == 0) {
	dump_mkdir("fields");
	dump_mkdir("hydro");
	dump_mkdir("rundata");
	dump_mkdir("injectors");
	dump_mkdir("restart1");  // 1st backup
	dump_mkdir("restart2");  // 2nd backup
	dump_mkdir("particle");
	
	//dump_grid("rundata/grid");
	//dump_materials("rundata/materials");
	//dump_species("rundata/species");
	global_header(g, np, "global", diag_.outputParams);
      }

      // Normal rundata energies dump
      if(should_dump(energies)) {
      	dump_energies("rundata/energies", step != 0, mprts, mflds, interpolator);
      }
      
      // Field data output

      if(should_dump(fields)) field_dump(mflds, diag_.fdParams);
      
      // Species moment output
      
      if(should_dump(ehydro)) hydro_dump(mprts, interpolator, mflds_hydro, "electron", diag_.hedParams);
      if(should_dump(Hhydro)) hydro_dump(mprts, interpolator, mflds_hydro, "ion", diag_.hHdParams);

#if 0
      if(step && !(step % diag->restart_interval)) {
	if(!diag->rtoggle) {
	  diag->rtoggle = 1;
	  //checkpt("restart1/restart", 0);
	}
	else {
	  diag->rtoggle = 0;
	  //checkpt("restart2/restart", 0);
	} // if
      } // if
#endif
      
      // Dump particle data

#if 0
      char subdir[36];
      if ( should_dump(eparticle) && step !=0
	   && step > 56*(diag->fields_interval)  ) {
	// if ( should_dump(eparticle) && step !=0 ) {
	sprintf(subdir,"particle/T.%lld",step);
	dump_mkdir(subdir);
	sprintf(subdir,"particle/T.%lld/eparticle",step);
	simulation->dump_particles("electron", subdir);
      }

      if ( should_dump(Hparticle) && step !=0
	   && step > 56*(diag->fields_interval)  ) {
	sprintf(subdir,"particle/T.%lld/Hparticle",step);
	simulation->dump_particles("ion", subdir);
      }
#endif
    } TOC(user_diagnostics, 1);
  }

#undef should_dump

  void create_field_list(char * strlist, DumpParameters & dumpParams)
  {
    strcpy(strlist, "");
    for(int i = 0, pass = 0; i < total_field_groups; i++) {
      if(dumpParams.output_vars.bitset(field_indeces[i])) {
	if(i>0 && pass) strcat(strlist, ", ");
	else pass = 1;
	strcat(strlist, fieldInfo[i].name);
      }
    }
  }
  
  void create_hydro_list(char * strlist, DumpParameters & dumpParams)
  {
    strcpy(strlist, "");
    for(size_t i(0), pass(0); i<total_hydro_groups; i++) {
      if(dumpParams.output_vars.bitset(hydro_indeces[i])) {
	if(i>0 && pass) strcat(strlist, ", ");
	else pass = 1;
	strcat(strlist, hydroInfo[i].name);
      }
    }
  }

  static int dump_mkdir(const char * dname)
  {
    return FileUtils::makeDirectory(dname);
  }

  static void print_hashed_comment( FileIO & fileIO, const char * comment)
  {
    fileIO.print("################################################################################\n");
    fileIO.print("# %s\n", comment);
    fileIO.print("################################################################################\n");
  }

  // ----------------------------------------------------------------------
  // global_header

  void global_header(const Grid *grid, const int *np, const char * base, std::vector<DumpParameters*> dumpParams)
  {
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) return;

    // Open the file for output
    char filename[256];
    sprintf(filename, "%s.vpc", base);
    
    FileIO fileIO;
    FileIOStatus status;
    
    status = fileIO.open(filename, io_write);
    if(status == fail) LOG_ERROR("Failed opening file: %s", filename);
    
    print_hashed_comment(fileIO, "Header version information");
    fileIO.print("VPIC_HEADER_VERSION 1.0.0\n\n");
    
    print_hashed_comment(fileIO,
			 "Header size for data file headers in bytes");
    fileIO.print("DATA_HEADER_SIZE 123\n\n");
    
    // Global grid inforation
    print_hashed_comment(fileIO, "Time step increment");
    fileIO.print("GRID_DELTA_T %f\n\n", grid->dt);
    
    print_hashed_comment(fileIO, "GRID_CVAC");
    fileIO.print("GRID_CVAC %f\n\n", grid->cvac);
    
    print_hashed_comment(fileIO, "GRID_EPS0");
    fileIO.print("GRID_EPS0 %f\n\n", grid->eps0);
    
    print_hashed_comment(fileIO, "Grid extents in the x-dimension");
    fileIO.print("GRID_EXTENTS_X %f %f\n\n", grid->x0, grid->x1);
    
    print_hashed_comment(fileIO, "Grid extents in the y-dimension");
    fileIO.print("GRID_EXTENTS_Y %f %f\n\n", grid->y0, grid->y1);
    
    print_hashed_comment(fileIO, "Grid extents in the z-dimension");
    fileIO.print("GRID_EXTENTS_Z %f %f\n\n", grid->z0, grid->z1);
    
    print_hashed_comment(fileIO, "Spatial step increment in x-dimension");
    fileIO.print("GRID_DELTA_X %f\n\n", grid->dx);
    
    print_hashed_comment(fileIO, "Spatial step increment in y-dimension");
    fileIO.print("GRID_DELTA_Y %f\n\n", grid->dy);
    
    print_hashed_comment(fileIO, "Spatial step increment in z-dimension");
    fileIO.print("GRID_DELTA_Z %f\n\n", grid->dz);
    
    print_hashed_comment(fileIO, "Domain partitions in x-dimension");
    fileIO.print("GRID_TOPOLOGY_X %d\n\n", np[0]);
    
    print_hashed_comment(fileIO, "Domain partitions in y-dimension");
    fileIO.print("GRID_TOPOLOGY_Y %d\n\n", np[1]);
    
    print_hashed_comment(fileIO, "Domain partitions in z-dimension");
    fileIO.print("GRID_TOPOLOGY_Z %d\n\n", np[2]);
    
    // Global data information
    assert(dumpParams.size() >= 2);
    
    print_hashed_comment(fileIO, "Field data information");
    fileIO.print("FIELD_DATA_DIRECTORY %s\n", dumpParams[0]->baseDir);
    fileIO.print("FIELD_DATA_BASE_FILENAME %s\n",
		 dumpParams[0]->baseFileName);
    
    // Create a variable list of field values to output.
    size_t numvars = std::min(dumpParams[0]->output_vars.bitsum(field_indeces,
								total_field_groups),
			      total_field_groups);
    size_t * varlist = new size_t[numvars];
    for(size_t v(0), c(0); v<total_field_groups; v++)
      if(dumpParams[0]->output_vars.bitset(field_indeces[v]))
	varlist[c++] = v;
    
    // output variable list
    fileIO.print("FIELD_DATA_VARIABLES %d\n", numvars);
    
    for(size_t v(0); v<numvars; v++)
      fileIO.print("\"%s\" %s %s %s %d\n", fieldInfo[varlist[v]].name,
		   fieldInfo[varlist[v]].degree, fieldInfo[varlist[v]].elements,
		   fieldInfo[varlist[v]].type, fieldInfo[varlist[v]].size);
    
    fileIO.print("\n");
    
    delete[] varlist;
    varlist = NULL;
    
    // Create a variable list for each species to output
    print_hashed_comment(fileIO, "Number of species with output data");
    fileIO.print("NUM_OUTPUT_SPECIES %d\n\n", dumpParams.size()-1);
    char species_comment[128];
    for(size_t i(1); i<dumpParams.size(); i++) {
      numvars = std::min(dumpParams[i]->output_vars.bitsum(hydro_indeces,
							   total_hydro_groups),
			 total_hydro_groups);
      
      sprintf(species_comment, "Species(%d) data information", (int)i);
      print_hashed_comment(fileIO, species_comment);
      fileIO.print("SPECIES_DATA_DIRECTORY %s\n",
		   dumpParams[i]->baseDir);
      fileIO.print("SPECIES_DATA_BASE_FILENAME %s\n",
		   dumpParams[i]->baseFileName);
      
      fileIO.print("HYDRO_DATA_VARIABLES %d\n", numvars);
      
      varlist = new size_t[numvars];
      for(size_t v(0), c(0); v<total_hydro_groups; v++)
	if(dumpParams[i]->output_vars.bitset(hydro_indeces[v]))
	  varlist[c++] = v;
      
      for(size_t v(0); v<numvars; v++)
	fileIO.print("\"%s\" %s %s %s %d\n", hydroInfo[varlist[v]].name,
		     hydroInfo[varlist[v]].degree, hydroInfo[varlist[v]].elements,
		     hydroInfo[varlist[v]].type, hydroInfo[varlist[v]].size);
      
      
      delete[] varlist;
      varlist = NULL;
      
      if(i<dumpParams.size()-1) fileIO.print("\n");
    }
    
    
    if( fileIO.close() ) LOG_ERROR("File close failed on global header!!!");
  }
  
  // ----------------------------------------------------------------------
  // dump_energies

  void dump_energies(const char *fname, int append,
		     Mparticles& mprts, MfieldsState& mflds, MfieldsInterpolator& interpolator)
  {
    double en_f[6], en_p;
    const Grid* g = mflds.vgrid();
    FileIO fileIO;
    FileIOStatus status(fail);

    int rank = psc_world_rank;
    
    if (!fname) LOG_ERROR("Invalid file name");
 
    if (rank==0) {
      status = fileIO.open(fname, append ? io_append : io_write);
      if (status==fail) LOG_ERROR("Could not open \"%s\".", fname);
      else {
	if (append==0) {
	  fileIO.print("%% Layout\n%% step ex ey ez bx by bz" );
	  for (auto& sp : mprts) {
	    fileIO.print(" \"%s\"", sp.name);
	  }
	  fileIO.print("\n");
	  fileIO.print("%% timestep = %e\n", g->dt);
	}
	fileIO.print( "%li", g->step );
      }
    }

    DiagOps::energy_f(mflds, en_f);
    if (rank==0 && status!=fail )
      fileIO.print( " %e %e %e %e %e %e",
		    en_f[0], en_f[1], en_f[2],
		    en_f[3], en_f[4], en_f[5] );
 
    for (auto& sp : mprts) {
      en_p = ParticlesOps::energy_p(sp, interpolator);
      if (rank==0 && status != fail) fileIO.print(" %e", en_p);
    }
 
    if (rank==0 && status != fail) {
      fileIO.print("\n");
      if (fileIO.close()) LOG_ERROR("File close failed on dump energies!!!");
    }
  }

  // ----------------------------------------------------------------------
  // field_dump
  
  void field_dump(MfieldsState& mflds, DumpParameters& dumpParams)
  {
    const Grid* grid = mflds.vgrid();
    auto& fa = mflds.getPatch(0);
    Field3D<typename MfieldsState::Patch> F(fa);
    int64_t step = grid->step;
    int rank = psc_world_rank;
    int nproc = psc_world_size;
    
    // Create directory for this time step
    char timeDir[256];
    sprintf(timeDir, "%s/T.%ld", dumpParams.baseDir, (long) step);
    dump_mkdir(timeDir);
  
    // Open the file for output
    char filename[256];
    sprintf(filename, "%s/T.%ld/%s.%ld.%d", dumpParams.baseDir, (long) step,
	    dumpParams.baseFileName, (long)step, rank);
  
    FileIO fileIO;
    FileIOStatus status;
  
    status = fileIO.open(filename, io_write);
    if( status==fail ) LOG_ERROR("Failed opening file: %s", filename);
  
    // convenience
    const size_t istride(dumpParams.stride_x);
    const size_t jstride(dumpParams.stride_y);
    const size_t kstride(dumpParams.stride_z);
  
    // Check stride values.
    if(remainder(grid->nx, istride) != 0)
      LOG_ERROR("x stride must be an integer factor of nx");
    if(remainder(grid->ny, jstride) != 0)
      LOG_ERROR("y stride must be an integer factor of ny");
    if(remainder(grid->nz, kstride) != 0)
      LOG_ERROR("z stride must be an integer factor of nz");
  
    int dim[3];
  
    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    int nxout = (grid->nx)/istride;
    int nyout = (grid->ny)/jstride;
    int nzout = (grid->nz)/kstride;
    float dxout = (grid->dx)*istride;
    float dyout = (grid->dy)*jstride;
    float dzout = (grid->dz)*kstride;

    assert(dumpParams.format == band);
    
    _WRITE_HEADER_V0(dump_type::field_dump, -1, 0, fileIO);
    
    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;
    
    WRITE_ARRAY_HEADER((&F(0,0,0)), 3, dim, fileIO);
    
    // Create a variable list of field values to output.
    size_t numvars = std::min(dumpParams.output_vars.bitsum(),
			      total_field_variables);
    size_t * varlist = new size_t[numvars];
    
    for(size_t i(0), c(0); i<total_field_variables; i++)
      if(dumpParams.output_vars.bitset(i)) varlist[c++] = i;
    
    if (numvars > 20) numvars = 20; // FIXME!!! materialid are 16 bit

    // more efficient for standard case
    if(istride == 1 && jstride == 1 && kstride == 1)
      for(size_t v(0); v<numvars; v++) {
	for(size_t k(0); k<nzout+2; k++) {
	  for(size_t j(0); j<nyout+2; j++) {
	    for(size_t i(0); i<nxout+2; i++) {
	      const uint32_t * fref = reinterpret_cast<uint32_t *>(&F(i,j,k));
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
	      const uint32_t * fref = reinterpret_cast<uint32_t *>(&F(ioff,joff,koff));
	      fileIO.write(&fref[varlist[v]], 1);
	    }
	  }
	}
      }
    
    delete[] varlist;
    
    if( fileIO.close() ) LOG_ERROR("File close failed on field dump!!!");
  }

  // ----------------------------------------------------------------------
  // hydro_dump
  
  void hydro_dump(Mparticles& mprts, MfieldsInterpolator& interpolator, MfieldsHydro& mflds_hydro,
		  const char * speciesname, DumpParameters& dumpParams )
  {
    Field3D<typename MfieldsHydro::Patch> H{mflds_hydro.getPatch(0)};
    const Grid* grid = mflds_hydro.vgrid();
    int64_t step = grid->step;
    int rank = psc_world_rank;
    int nproc = psc_world_size;

    // Create directory for this time step
    char timeDir[256];
    sprintf(timeDir, "%s/T.%ld", dumpParams.baseDir, (long)step);
    dump_mkdir(timeDir);
  
    // Open the file for output
    char filename[256];
    sprintf( filename, "%s/T.%ld/%s.%ld.%d", dumpParams.baseDir, (long)step,
	     dumpParams.baseFileName, (long)step, rank );
  
    FileIO fileIO;
    FileIOStatus status;

    status = fileIO.open(filename, io_write);
    if(status == fail) LOG_ERROR("Failed opening file: %s", filename);

    auto& sp = *std::find_if(mprts.begin(), mprts.end(),
			     [&](const typename Mparticles::Species &sp) { return strcmp(sp.name, speciesname) == 0; });

    HydroArrayOps::clear(mflds_hydro);
    ParticlesOps::accumulate_hydro_p(mflds_hydro, sp, interpolator);
    HydroArrayOps::synchronize(mflds_hydro);
  
    // convenience
    const size_t istride(dumpParams.stride_x);
    const size_t jstride(dumpParams.stride_y);
    const size_t kstride(dumpParams.stride_z);
  
    // Check stride values.
    if(remainder(grid->nx, istride) != 0)
      LOG_ERROR("x stride must be an integer factor of nx");
    if(remainder(grid->ny, jstride) != 0)
      LOG_ERROR("y stride must be an integer factor of ny");
    if(remainder(grid->nz, kstride) != 0)
      LOG_ERROR("z stride must be an integer factor of nz");
  
    int dim[3];
  
    /* IMPORTANT: these values are written in WRITE_HEADER_V0 */
    int nxout = (grid->nx)/istride;
    int nyout = (grid->ny)/jstride;
    int nzout = (grid->nz)/kstride;
    float dxout = (grid->dx)*istride;
    float dyout = (grid->dy)*jstride;
    float dzout = (grid->dz)*kstride;
  
    assert(dumpParams.format == band);
    
    _WRITE_HEADER_V0(dump_type::hydro_dump, sp.id, sp.q/sp.m, fileIO);
    
    dim[0] = nxout+2;
    dim[1] = nyout+2;
    dim[2] = nzout+2;
    
    WRITE_ARRAY_HEADER((&H(0,0,0)), 3, dim, fileIO);
    
    /*
     * Create a variable list of hydro values to output.
     */
    size_t numvars = std::min(dumpParams.output_vars.bitsum(),
			      total_hydro_variables);
    size_t * varlist = new size_t[numvars];
    for(size_t i(0), c(0); i<total_hydro_variables; i++)
      if( dumpParams.output_vars.bitset(i) ) varlist[c++] = i;
    
    // More efficient for standard case
    if(istride == 1 && jstride == 1 && kstride == 1)
      
      for(size_t v(0); v<numvars; v++)
	for(size_t k(0); k<nzout+2; k++)
	  for(size_t j(0); j<nyout+2; j++)
	    for(size_t i(0); i<nxout+2; i++) {
	      const uint32_t * href = reinterpret_cast<uint32_t *>(&H(0,0,0));
	      fileIO.write(&href[varlist[v]], 1);
	    }
    
    else
      
      for(size_t v(0); v<numvars; v++)
	for(size_t k(0); k<nzout+2; k++) { const size_t koff = (k == 0) ? 0 : (k == nzout+1) ? grid->nz+1 : k*kstride-1;
	  for(size_t j(0); j<nyout+2; j++) { const size_t joff = (j == 0) ? 0 : (j == nyout+1) ? grid->ny+1 : j*jstride-1;
	    for(size_t i(0); i<nxout+2; i++) { const size_t ioff = (i == 0) ? 0 : (i == nxout+1) ? grid->nx+1 : i*istride-1;
	      const uint32_t * href = reinterpret_cast<uint32_t *>(&H(ioff,joff,koff));
	      fileIO.write(&href[varlist[v]], 1);
	    }
	  }
	}
    
    delete[] varlist;
    
    if( fileIO.close() ) LOG_ERROR("File close failed on hydro dump!!!");
  }

private:
  VpicDiag diag_;
};


#endif
