
#ifndef PSC_VPIC_BITS_H
#define PSC_VPIC_BITS_H

#include <mpi.h>

#include <cstdint>

extern MPI_Comm psc_comm_world;
extern int psc_world_rank;
extern int psc_world_size;

// FIXME, this is chosen globally and to match vpic, but could be
// made dependent on actual MaterialList implementation
enum { MaterialIdMax = 32768 };
typedef int16_t MaterialId;

// FIXME, this is chosen globally and to match vpic for now
typedef int32_t SpeciesId; // Must be 32-bit wide for particle_injector_t

// FIXME, get rid of this BOUNDARY macro
#define BOUNDARY(i,j,k) (13+(i)+3*(j)+9*(k)) /* FORTRAN -1:1,-1:1,-1:1 */

// FIXME, get rid of
#define VOXEL(x,y,z, nx,ny,nz) ((x) + ((nx)+2)*((y) + ((ny)+2)*(z)))

// FIXME, get rid of
#define BEGIN_PRIMITIVE do
#define END_PRIMITIVE   while(0)

// FIXME
#ifndef ALIGNED
#define ALIGNED(a)
#endif

// FIXME
#ifndef RESTRICT
#define RESTRICT __restrict
#endif 

// FIXME
#define DECLARE_ALIGNED_ARRAY(type,align,name,count)    \
  type name[(count)] __attribute__ ((aligned (align)))

// FIXME
#define POW2_CEIL(u,a) ( ((u)+(a)-1) & (~((a)-1)) )

// FIXME
#ifndef LIKELY
#define LIKELY(_c)   __builtin_expect((_c),1)
#endif
#ifndef UNLIKELY
#define UNLIKELY(_c) __builtin_expect((_c),0)
#endif

#define mprintf(fmt...) do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); { printf("[%d] ", __rank); printf(fmt); } } while(0)
#define MHERE do { int __rank; MPI_Comm_rank(MPI_COMM_WORLD, &__rank); printf("[%d] HERE: in %s() at %s:%d\n", __rank, __FUNCTION__, __FILE__, __LINE__); } while(0)

// FIXME, do some proper logging eventually...

#define LOG_ERROR(fmt...) do {						\
    mprintf("ERROR at %s:%d (%s): ", __FILE__, __LINE__, __func__); printf(fmt); \
    abort();								\
  } while (0)
    
#define LOG_WARN(fmt...) do {						\
    mprintf("WARNING at %s:%d (%s): ", __FILE__, __LINE__, __func__); printf(fmt); \
  } while (0)
    
#define LOG_INFO(fmt...) do {						\
    mprintf("INFO at %s:%d (%s): ", __FILE__, __LINE__, __func__); printf(fmt); \
  } while (0)
    

#endif
