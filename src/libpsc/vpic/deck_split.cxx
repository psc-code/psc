
#include "vpic_iface.h"

#include "vpic_init.h"

#include <mrc_common.h>

#include <cassert>

#undef sim_log
#define sim_log(x) do {                                \
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);    \
    if( rank==0 ) {				       \
      std::cerr << "SIM_LOG: " << x << std::endl;      \
      std::cerr.flush();                               \
    }                                                  \
  } while(0)

// ----------------------------------------------------------------------
// vpic_simulation_setup_diagnostics

void vpic_simulation_setup_diagnostics(vpic_simulation *simulation, globals_diag *diag,
				       species_t *electron, species_t *ion)
{
  /*--------------------------------------------------------------------------
   * New dump definition
   *------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------
   * Set data output format
   *
   * This option allows the user to specify the data format for an output
   * dump.  Legal settings are 'band' and 'band_interleave'.  Band-interleave
   * format is the native storage format for data in VPIC.  For field data,
   * this looks something like:
   *
   *   ex0 ey0 ez0 div_e_err0 cbx0 ... ex1 ey1 ez1 div_e_err1 cbx1 ...
   *
   * Banded data format stores all data of a particular state variable as a
   * contiguous array, and is easier for ParaView to process efficiently.
   * Banded data looks like:
   *
   *   ex0 ex1 ex2 ... exN ey0 ey1 ey2 ...
   *
   *------------------------------------------------------------------------*/

  diag->fdParams.format = band;
  sim_log( "Fields output format = band" );

  diag->hedParams.format = band;
  sim_log( "Electron species output format = band" );

  diag->hHdParams.format = band;
  sim_log( "Ion species output format = band" );

  /*--------------------------------------------------------------------------
   * Set stride
   *
   * This option allows data down-sampling at output.  Data are down-sampled
   * in each dimension by the stride specified for that dimension.  For
   * example, to down-sample the x-dimension of the field data by a factor
   * of 2, i.e., half as many data will be output, select:
   *
   *   global->fdParams.stride_x = 2;
   *
   * The following 2-D example shows down-sampling of a 7x7 grid (nx = 7,
   * ny = 7.  With ghost-cell padding the actual extents of the grid are 9x9.
   * Setting the strides in x and y to equal 2 results in an output grid of
   * nx = 4, ny = 4, with actual extents 6x6.
   *
   * G G G G G G G G G
   * G X X X X X X X G
   * G X X X X X X X G         G G G G G G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G   ==>   G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G X X X X G
   * G X X X X X X X G         G G G G G G
   * G G G G G G G G G
   *
   * Note that grid extents in each dimension must be evenly divisible by
   * the stride for that dimension:
   *
   *   nx = 150;
   *   global->fdParams.stride_x = 10; // legal -> 150/10 = 15
   *
   *   global->fdParams.stride_x = 8; // illegal!!! -> 150/8 = 18.75
   *------------------------------------------------------------------------*/

  // relative path to fields data from global header
  sprintf(diag->fdParams.baseDir, "fields");

  // base file name for fields output
  sprintf(diag->fdParams.baseFileName, "fields");

  diag->fdParams.stride_x = 1;
  diag->fdParams.stride_y = 1;
  diag->fdParams.stride_z = 1;

  // add field parameters to list
  diag->outputParams.push_back(&diag->fdParams);

  sim_log( "Fields x-stride " << diag->fdParams.stride_x );
  sim_log( "Fields y-stride " << diag->fdParams.stride_y );
  sim_log( "Fields z-stride " << diag->fdParams.stride_z );

  // relative path to electron species data from global header
  sprintf(diag->hedParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(diag->hedParams.baseFileName, "ehydro");

  diag->hedParams.stride_x = 1;
  diag->hedParams.stride_y = 1;
  diag->hedParams.stride_z = 1;

  // add electron species parameters to list
  diag->outputParams.push_back(&diag->hedParams);

  sim_log( "Electron species x-stride " << diag->hedParams.stride_x );
  sim_log( "Electron species y-stride " << diag->hedParams.stride_y );
  sim_log( "Electron species z-stride " << diag->hedParams.stride_z );

  // relative path to electron species data from global header
  sprintf(diag->hHdParams.baseDir, "hydro");

  // base file name for fields output
  sprintf(diag->hHdParams.baseFileName, "Hhydro");

  diag->hHdParams.stride_x = 1;
  diag->hHdParams.stride_y = 1;
  diag->hHdParams.stride_z = 1;

  sim_log( "Ion species x-stride " << diag->hHdParams.stride_x );
  sim_log( "Ion species y-stride " << diag->hHdParams.stride_y );
  sim_log( "Ion species z-stride " << diag->hHdParams.stride_z );

  // add electron species parameters to list
  diag->outputParams.push_back(&diag->hHdParams);

  /*--------------------------------------------------------------------------
   * Set output fields
   *
   * It is now possible to select which state-variables are output on a
   * per-dump basis.  Variables are selected by passing an or-list of
   * state-variables by name.  For example, to only output the x-component
   * of the electric field and the y-component of the magnetic field, the
   * user would call output_variables like:
   *
   *   global->fdParams.output_variables( ex | cby );
   *
   * NOTE: OUTPUT VARIABLES ARE ONLY USED FOR THE BANDED FORMAT.  IF THE
   * FORMAT IS BAND-INTERLEAVE, ALL VARIABLES ARE OUTPUT AND CALLS TO
   * 'output_variables' WILL HAVE NO EFFECT.
   *
   * ALSO: DEFAULT OUTPUT IS NONE!  THIS IS DUE TO THE WAY THAT VPIC
   * HANDLES GLOBAL VARIABLES IN THE INPUT DECK AND IS UNAVOIDABLE.
   *
   * For convenience, the output variable 'all' is defined:
   *
   *   global->fdParams.output_variables( all );
   *------------------------------------------------------------------------*/
  /* CUT AND PASTE AS A STARTING POINT
   * REMEMBER TO ADD APPROPRIATE GLOBAL DUMPPARAMETERS VARIABLE

   output_variables( all );

   output_variables( electric | div_e_err | magnetic | div_b_err |
                     tca      | rhob      | current  | rhof |
                     emat     | nmat      | fmat     | cmat );

   output_variables( current_density  | charge_density |
                     momentum_density | ke_density     | stress_tensor );
   */

  //global->fdParams.output_variables( electric | magnetic );
  //global->hedParams.output_variables( current_density | charge_density
  //                                    | stress_tensor );
  //global->hHdParams.output_variables( current_density | charge_density );
  //                                    | stress_tensor );

  diag->fdParams.output_variables( all );
  diag->hedParams.output_variables( all );
  diag->hHdParams.output_variables( all );

  /*--------------------------------------------------------------------------
   * Convenience functions for simlog output
   *------------------------------------------------------------------------*/

  char varlist[512];
  simulation->create_field_list(varlist, diag->fdParams);

  sim_log( "Fields variable list: " << varlist );

  simulation->create_hydro_list(varlist, diag->hedParams);

  sim_log( "Electron species variable list: " << varlist );

  simulation->create_hydro_list(varlist, diag->hHdParams);

  sim_log( "Ion species variable list: " << varlist );
}

#define should_dump(x)                                                  \
  (diag->x##_interval>0 && remainder(step, diag->x##_interval) == 0)

void vpic_simulation_diagnostics(vpic_simulation *simulation, globals_diag *diag)
{
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

  /*--------------------------------------------------------------------------
   * Normal rundata dump
   *------------------------------------------------------------------------*/
  if(step==0) {
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
    simulation->global_header("global", diag->outputParams);
  } // if

  /*--------------------------------------------------------------------------
   * Normal rundata energies dump
   *------------------------------------------------------------------------*/
  if(should_dump(energies)) {
    simulation->dump_energies("rundata/energies", step == 0 ? 0 : 1);
  } // if

  /*--------------------------------------------------------------------------
   * Field data output
   *------------------------------------------------------------------------*/

  if(step == -1 || should_dump(fields)) simulation->field_dump(diag->fdParams);

  /*--------------------------------------------------------------------------
   * Electron species output
   *------------------------------------------------------------------------*/

  if(should_dump(ehydro)) simulation->hydro_dump("electron", diag->hedParams);

  /*--------------------------------------------------------------------------
   * Ion species output
   *------------------------------------------------------------------------*/

  if(should_dump(Hhydro)) simulation->hydro_dump("ion", diag->hHdParams);

  /*--------------------------------------------------------------------------
   * Restart dump
   *------------------------------------------------------------------------*/

  if(step && !(step%diag->restart_interval)) {
    if(!diag->rtoggle) {
      diag->rtoggle = 1;
      checkpt("restart1/restart", 0);
    }
    else {
      diag->rtoggle = 0;
      checkpt("restart2/restart", 0);
    } // if
  } // if

  // Dump particle data

  char subdir[36];
  if ( should_dump(eparticle) && step !=0
       && step > 56*(diag->fields_interval)  ) {
    // if ( should_dump(eparticle) && step !=0 ) {
    sprintf(subdir,"particle/T.%lld",step);
    simulation->dump_mkdir(subdir);
    sprintf(subdir,"particle/T.%lld/eparticle",step);
    simulation->dump_particles("electron", subdir);
  }

  if ( should_dump(Hparticle) && step !=0
       && step > 56*(diag->fields_interval)  ) {
    sprintf(subdir,"particle/T.%lld/Hparticle",step);
    simulation->dump_particles("ion", subdir);
  }

}

