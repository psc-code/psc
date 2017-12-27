
#ifndef PSC_VPIC_BITS_H
#define PSC_VPIC_BITS_H

#include <mpi.h>

extern MPI_Comm psc_comm_world;
extern int psc_world_rank;
extern int psc_world_size;

#include "vpic.h" // FIXME

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)
#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

// FIXME, do some proper logging eventually...

#define LOG_ERROR(fmt...) do {						\
    mprintf("ERROR at %s:%d (%s): ", __FILE__, __LINE__, __func__); printf(fmt); \
    abort();								\
  } while (0)
    
#define LOG_WARN(fmt...) do {						\
    mprintf("WARNING at %s:%d (%s): ", __FILE__, __LINE__, __func__); printf(fmt); \
    abort();								\
  } while (0)
    
#define LOG_INFO(fmt...) do {						\
    mprintf("INFO at %s:%d (%s): ", __FILE__, __LINE__, __func__); printf(fmt); \
    abort();								\
  } while (0)
    

#endif
