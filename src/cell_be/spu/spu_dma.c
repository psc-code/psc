#include "../psc_cbe_common.h"
#ifdef __SPU__
#include <spu_mfcio.h>
#else
#include "spu_mfcio_c.h"
#endif

#include <assert.h>

enum {
  ST_NULL, 
  ST_LOADING, 
  ST_LOADED, 
  ST_STORING, 
  ST_STORED,
};

struct bl_info {
  int state; 
  int tagid; 
};

// Is this a hardware limit?
#define MAX_TAGID (32)


struct mfc_dma_el {
  unsigned int size;
  unsigned long eal;
};

static unsigned int tagmask, tagmask_storing, tagmask_loading;

// returns the first available tag, and sets
// that bit to 1 in the tagmask 
static int
get_tagid(void)
{
  for (int i = 0; i < MAX_TAGID; i++) {
    if (!(tagmask & (1 << i))) {
      tagmask |= (1 << i);
      return i;
    }
  }
  return -1;
}

// sets the tagid-th bit of tagmask 
// to 0
static void
put_tagid(int tagid)
{
  assert(tagid >= 0 && tagid < MAX_TAGID);
  assert(tagmask & (1 << tagid));
  tagmask &= ~(1 << tagid);
}

	

static void
wait_tagid(int tagid)
{
  assert(tagid >= 0 && tagid < MAX_TAGID);
  assert(tagmask & (1 << tagid));
  // Assign what tagid to look for 
  mfc_write_tag_mask(1 << tagid);
  // Wait until the previously indicated tagid
  // is done transfering. 
  unsigned int mask = mfc_read_tag_status_any();

  // the openggcm implementation has 
  // some other stuff here. I'm not 
  // sure it's needed for psc. 

}

void
spu_dma_get(volatile void *ls, unsigned long long ea, unsigned long size)
{
  // Check that we're on 16B boundaries, and 
  // the size of the struct we're bringing in is 
  // a multiple of 16B 
  //  fprintf(stderr, "size %d\n", size);
  assert(((unsigned long)ls & 15) == 0);
  assert((ea & 15) == 0);
  assert((size & 15) == 0);
  //fflush(stdout);
  //  fprintf(stderr,"dma_get %p %#llx %#lx\n", ls, ea, size);

  int tagid = get_tagid();
  assert(tagid >= 0);
  
  //  fprintf(stderr, " size %d \n", size);
  mfc_get(ls, ea, size, tagid, 0, 0);
  wait_tagid(tagid);
  put_tagid(tagid);
}
