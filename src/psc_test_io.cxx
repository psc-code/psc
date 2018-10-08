
#include "vec3.hxx"
#include "psc_test_io_xdmf.h"

#include <mrc_domain.h>

#include <string>

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

// ======================================================================
// PscTestIo

struct PscTestIo
{
  static void run()
  {
    mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

    // --- setup domain
#if 0
    Int3 gdims = { 400, 800, 2400}; // global number of grid points
    Int3 np = { 8, 16, 48 }; // division into patches
#else
    Int3 gdims = { 20, 20, 80}; // global number of grid points
    Int3 np = { 2, 2, 8 }; // division into patches
#endif
    
    struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
    mrc_domain_set_type(domain, "multi");
    mrc_domain_set_param_int3(domain, "m", gdims);
    mrc_domain_set_param_int(domain, "bcx", BC_PERIODIC);
    mrc_domain_set_param_int(domain, "bcy", BC_PERIODIC);
    mrc_domain_set_param_int(domain, "bcz", BC_PERIODIC);
    mrc_domain_set_param_int3(domain, "np", np);
    
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    
    mpi_printf(MPI_COMM_WORLD, "***** Testing output\n");

    mrc_fld* fld = mrc_domain_m3_create(domain);
    mrc_fld_set_name(fld, "e");
    mrc_fld_set_param_int(fld, "nr_ghosts", 0);
    mrc_fld_set_param_int(fld, "nr_comps", 2);
    mrc_fld_setup(fld);
    mrc_fld_set_comp_name(fld, 0, "ex");
    mrc_fld_set_comp_name(fld, 1, "ey");
    
    for (int p = 0; p < mrc_fld_nr_patches(fld); p++) {
      mrc_fld_patch *m3p = mrc_fld_patch_get(fld, p);
      mrc_fld_foreach(fld, i,j,k, 0,0) {
	MRC_M3(m3p, 0, i,j,k) = i;
	MRC_M3(m3p, 1, i,j,k) = j;
      } mrc_fld_foreach_end;
      mrc_fld_patch_put(fld);
    }
    
    xdmf xdmf[1] = {};
    xdmf->nr_writers = 2;
    xdmf_collective_setup(xdmf);
    xdmf_collective_write_m3(xdmf, "testpath", fld);
    xdmf_collective_destroy(xdmf);

    mrc_fld_destroy(fld);
    mrc_domain_destroy(domain);

    mpi_printf(MPI_COMM_WORLD, "***** Testing output done\n");
  }
};


// ======================================================================
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  PscTestIo::run();

  MPI_Finalize();
  return 0;
}
