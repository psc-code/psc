
#include "dim.hxx"
#include "grid.hxx"
#include "psc_fields_single.h"

#include "mrc_io_private.h"

#include <string>

static inline struct mrc_io_ops *
mrc_io_ops(struct mrc_io *io)
{
  return (struct mrc_io_ops *) io->obj.ops;
}


int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
	      int tag, MPI_Comm comm, MPI_Request *request)
{
  mprintf("MPI_Irecv buf %p count %d source %d tag %d request %p\n",
	  buf, count, source, tag, request);
  int ierr = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  return ierr;
}


int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest,
	      int tag, MPI_Comm comm, MPI_Request *request)
{
  mprintf("MPI_Isend buf %p count %d dest %d tag %d request %p\n",
	  buf, count, dest, tag, request);
  int ierr = PMPI_Irecv(buf, count, datatype, dest, tag, comm, request);
  return ierr;
}

int MPI_Waitall(int count, MPI_Request array_of_requests[],
		MPI_Status array_of_statuses[])
{
  mprintf("MPI_Waitall count %d req %p status %p\n", count, array_of_requests, array_of_statuses);
  int ierr = PMPI_Waitall(count, array_of_requests, array_of_statuses);
  return ierr;
}

static unsigned long
mrc_obj_uid(struct mrc_obj *obj)
{
  unsigned long uid = (unsigned long) obj;
  if (obj->comm != MPI_COMM_NULL) {
    MPI_Bcast(&uid, 1, MPI_LONG, 0, obj->comm);
  }
  return uid;
}

// ======================================================================
// PscTestIo

struct PscTestIo
{
  static void run()
  {
    mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

    // -- setup particle kinds
    Grid_t::Kinds kinds = {};
    
    // --- setup domain
    Grid_t::Real3 LL = { 400., 800., 400.*6 }; // domain size (in d_e)
#if 0
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 8, 16, 48 }; // division into patches
#else
    Int3 gdims = { 20, 20, 80}; // global number of grid points
    Int3 np = { 2, 2, 8 }; // division into patches
#endif
    
    auto grid_domain = Grid_t::Domain{gdims, LL, -.5 * LL, np};
    
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};

    // --- generic setup
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 5;
    
    double dt = .99;
    auto norm = Grid_t::Normalization{norm_params};
    auto grid = Grid_t{grid_domain, grid_bc, kinds, norm, dt};

    mpi_printf(MPI_COMM_WORLD, "***** Testing output\n");

    auto io = mrc_io_create(MPI_COMM_WORLD);
    mrc_io_set_param_string(io, "basename", "pfd");
    mrc_io_set_param_string(io, "outdir", ".");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
    mrc_io_view(io);

    mrc_io_open(io, "w", 0, 0.);
    
    // save some basic info about the run in the output file
    struct mrc_obj *obj = mrc_obj_create(mrc_io_comm(io));
    mrc_obj_set_name(obj, "psc");
    mrc_obj_dict_add_int(obj, "timestep", 0);
    mrc_obj_dict_add_float(obj, "time", 0.);
    mrc_obj_dict_add_float(obj, "cc", 100.);
    mrc_obj_dict_add_float(obj, "dt", .99);
    mrc_obj_write(obj, io);
    mrc_obj_destroy(obj);
    
    mrc_fld* fld = grid.mrc_domain().m3_create();
    mrc_fld_set_name(fld, "e");
    mrc_fld_set_param_int(fld, "nr_ghosts", 0);
    mrc_fld_set_param_int(fld, "nr_comps", 2);
    mrc_fld_setup(fld);
    mrc_fld_set_comp_name(fld, 0, "ex");
    mrc_fld_set_comp_name(fld, 1, "ey");
    
    for (int p = 0; p < grid.n_patches(); p++) {
      mrc_fld_patch *m3p = mrc_fld_patch_get(fld, p);
      mrc_fld_foreach(fld, i,j,k, 0,0) {
	MRC_M3(m3p, 0, i,j,k) = i;
	MRC_M3(m3p, 1, i,j,k) = j;
      } mrc_fld_foreach_end;
      mrc_fld_patch_put(fld);
    }
    
    char path[strlen(obj->name) + 20];
    sprintf(path, "%s-uid-%#lx", obj->name, mrc_obj_uid(obj));
    if (mrc_io_add_obj(io, obj, path) == 1) {
      MHERE;
    }
    struct mrc_io_ops *ops = mrc_io_ops(io);
    ops->write_m3(io, path, fld);

    mrc_fld_destroy(fld);
    mrc_io_close(io);
    mrc_io_destroy(io);

    mpi_printf(MPI_COMM_WORLD, "***** Testing output done\n");
  }
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  PscTestIo::run();

  libmrc_params_finalize();
  MPI_Finalize();
  return 0;
}
