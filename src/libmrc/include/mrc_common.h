
#ifndef MRC_COMMON_H
#define MRC_COMMON_H

#include <config.h>

enum {
  BC_NONE,
  BC_PERIODIC,
};

#ifndef NDEBUG

#define assert_collective(comm) do {		\
    MPI_Barrier(comm);			\
  } while (0)

#else

#define assert_collective(comm) do {} while (0)

#endif

#define mpi_printf(comm, args...) do { int __rank; MPI_Comm_rank(comm, &__rank); if (__rank == 0) { printf(args); } } while(0)

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

#ifndef __unused
#define __unused __attribute__((__unused__))
#endif

#endif
