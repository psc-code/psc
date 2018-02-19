
#ifndef MRC_COMMON_H
#define MRC_COMMON_H

#include <mrc_config.h>
#include <mrc.h>

#include <mpi.h>

enum {
  BC_NONE,
  BC_PERIODIC,
};

#ifndef NDEBUG

#define assert_collective(comm) do {		\
    if (comm != MPI_COMM_NULL) {		\
      MPI_Barrier(comm);			\
    }						\
  } while (0)

#else

#define assert_collective(comm) do {} while (0)

#endif

extern int mrc_view_level;

void mrc_view_printf(MPI_Comm comm, const char *fmt, ...);

#define mpi_printf(comm, args...) do { int __rank; MPI_Comm_rank(comm, &__rank); if (__rank == 0) { printf(args); } } while(0)

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)

#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

#ifdef __GNUC__
#define	_mrc_unused	__attribute__((__unused__))
#else
#define	_mrc_unused	/* no attribute */
#endif

#ifndef __deprecated
  #ifdef __GNUC__
    #define __deprecated(func) (func) __attribute__ ((deprecated))
  #else
    #pragma message("WARNING: __deprecated not implemented for this compiler")
    #define __deprecated(func) (func)
  #endif
#endif

#endif
