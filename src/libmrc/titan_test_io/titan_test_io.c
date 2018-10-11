
#define BOUNDS_CHECK
#define BOUNDSCHECK

#include <mpi.h>

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)


int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // --- setup domain
#if 1
  int gdims[3] = { 400, 800, 2400}; // global number of grid points
  int np[3] = { 8, 16, 48 }; // division into patches
#else
  int gdims[3] = { 20, 20, 80}; // global number of grid points
  int np[3] = { 2, 2, 8 }; // division into patches
#endif

  int n_writers = 2;

  int ldims[3];
  for (int d = 0; d < 3; d++) {
    assert(gdims[d] % np[d] == 0);
    ldims[d] = gdims[d] / np[d];
  }

  int n_patches = np[0] * np[1] * np[2];
  assert(n_patches % size == 0);
  int n_patches_per_rank = n_patches / size;
  mprintf("n_patches %d per rank %d\n", n_patches, n_patches_per_rank);

  assert(gdims[2] % n_writers == 0);
  int writer_dims[3] = { gdims[0], gdims[1], gdims[2] / n_writers };
  int writer_offs[3] = {};
  if (rank < n_writers) {
    writer_offs[2] = rank * writer_dims[2];
  }

  size_t buf_size = n_patches_per_rank * ldims[0] * ldims[1] * ldims[2];
  float *send_buf = calloc(buf_size, sizeof(*send_buf));
  float *recv_buf = calloc(buf_size, sizeof(*recv_buf));
  MPI_Request *send_reqs = calloc(n_writers, sizeof(*send_reqs));
  MPI_Request *recv_reqs = NULL;
  int n_recv_reqs = 0;

  for (int m = 0; m < 2; m++) {
    MPI_Barrier(comm);
    if (rank == 0) printf("round %d\n", m);
    MPI_Barrier(comm);

    int patch = rank * n_patches_per_rank;
    int p = rank * n_patches_per_rank;
    int px = p % np[0]; p /= np[0];
    int py = p % np[1]; p /= np[1];
    int pz = p;

    int npz_writer = np[2] / n_writers;
    int writer = pz / npz_writer;
    mprintf("patch %d %d:%d:%d writer %d npz_writer %d\n", patch, px, py, pz, writer, npz_writer);

    if (rank < n_writers) {
      n_recv_reqs = size / n_writers;
      recv_reqs = calloc(n_recv_reqs, sizeof(*recv_reqs));
      for (int source = rank * n_recv_reqs; source < (rank + 1) * n_recv_reqs; source++) {
	mprintf("irecv src %d\n", source);
	MPI_Irecv(recv_buf, buf_size, MPI_FLOAT, source, 111, comm, &recv_reqs[source - rank * n_recv_reqs]);
      }
    } 
    for (int dest = 0; dest < n_writers; dest++) {
      if (dest == writer) {
	mprintf("isend dest %d\n", dest);
	MPI_Isend(send_buf, buf_size, MPI_FLOAT, dest, 111, comm, &send_reqs[dest]);
      } else {
	send_reqs[dest] = MPI_REQUEST_NULL;
      }
    }

    if (rank < n_writers) {
      MPI_Waitall(n_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE);
      free(recv_reqs);
    }
    
    MPI_Waitall(n_writers, send_reqs, MPI_STATUSES_IGNORE);
  }

  free(send_buf);
  free(send_reqs);
  free(recv_buf);
  
  MPI_Finalize();
  return 0;
}

