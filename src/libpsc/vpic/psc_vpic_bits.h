
#ifndef PSC_VPIC_BITS_H
#define PSC_VPIC_BITS_H

#include <mpi.h>

extern MPI_Comm psc_comm_world;
extern int psc_world_rank;
extern int psc_world_size;

#include "vpic.h" // FIXME

#undef ERROR
#define ERROR(args) do {						\
    mprintf( "Error at " _LOG_HDR ":\n\t" );				\
    mprintf args;							\
    mprintf( "\n" );							\
    exit(1);								\
  } while(0)

#undef WARNING
#define WARNING(args) do {						\
    mprintf( "Warning at " _LOG_HDR ":\n\t" );				\
    mprintf args;							\
    mprintf( "\n" );							\
  } while(0)


#endif
