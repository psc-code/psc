
#include "psc_test_io_xdmf.h"

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,
	      int tag, MPI_Comm comm, MPI_Request *request)
{
  mprintf("MPI_Irecv buf %p count %d source %d tag %d request %p\n",
	  buf, count, source, tag, request);
  int ierr = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
  return ierr;
}


int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
	      int tag, MPI_Comm comm, MPI_Request *request)
{
  mprintf("MPI_Isend buf %p count %d dest %d tag %d request %p\n",
	  buf, count, dest, tag, request);
  int ierr = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
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
    int gdims[3] = { 400, 800, 2400}; // global number of grid points
    int np[3] = { 8, 16, 48 }; // division into patches
#else
    int gdims[3] = { 20, 20, 80}; // global number of grid points
    int np[3] = { 2, 2, 8 }; // division into patches
#endif
    
    mpi_printf(MPI_COMM_WORLD, "***** Testing output\n");

    xdmf xdmf[1];
    xdmf_collective_setup(xdmf, 2, gdims, np);
    xdmf_collective_write_m3(xdmf);
    xdmf_collective_destroy(xdmf);

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
