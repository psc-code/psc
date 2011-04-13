
#include "spu_mfcio_c.h"
#include "../libspe2_c.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

// fake cell implementation...

static int tagmask;
static int __mask;

void
mfc_write_tag_mask(unsigned int mask)
{
  __mask = mask;
}

unsigned int
mfc_read_tag_status_any()
{
  unsigned int retval = tagmask & __mask;
  assert(retval);
  tagmask &= ~__mask;

  return retval;
}

struct mfc_dma_el { 
  unsigned int size; 
  unsigned long eal; 
};

void
mfc_get(volatile void *ls, unsigned long long ea, unsigned long size, 
	unsigned int tag, unsigned int tid, unsigned int rid)
{
  memcpy((void *) ls, (void *)(unsigned long) ea, size);
  tagmask |= (1 << tag);
}

void
mfc_getl(volatile void *ls, unsigned long long ea, void *lsa, 
	 unsigned long size, unsigned int tag, unsigned int tid, unsigned int rid)
{
  for (struct mfc_dma_el *p = lsa; p < (struct mfc_dma_el *)(lsa + size); p++) {
    assert((p->size & 15) == 0);
    memcpy((void *) ls, (void *) p->eal, p->size);
    ls += p->size;
  }
  tagmask |= (1 << tag);
}

extern float fldsp[];

void
mfc_putl(volatile void *ls, unsigned long long ea, void *lsa, 
	 unsigned long size, unsigned int tag, unsigned int tid, unsigned int rid)
{
  for (struct mfc_dma_el *p = lsa; p < (struct mfc_dma_el *)(lsa + size); p++) {
    if ((p->size & 15) == 0) {
      memcpy((void *) p->eal, (void *) ls, p->size);
      /* for (int i = 0; i < 4; i++) { */
      /* 	((float *) p->eal)[i] = i; */
      /* } */
      ls += p->size;
    } else if (p->size == 8 || p->size == 4) { // FIXME shouldn't occur (codegen side)
      memcpy((void *) p->eal, (void *) ((unsigned long)ls | (p->eal & 15)), p->size);
      ls += 16;
    } else {
      printf("p->size %d !!!\n", p->size);
      assert(0);
    }
  }
  tagmask |= (1 << tag);
}

// ----------------------------------------------------------------------
// mbox

#define MAX_SPES (8)

static pthread_mutex_t msg_mutex[MAX_SPES];
static pthread_cond_t msg_cond[MAX_SPES];
static int mbox_in[MAX_SPES];
static pthread_mutex_t mbox_out_mutex[MAX_SPES];
static int mbox_out[MAX_SPES];

// FIXME, instead we should put the above into TLS.
static pthread_key_t key; 

unsigned int
spu_read_in_mbox()
{
  int tid = (int) (unsigned long) pthread_getspecific(key);
  //  printf("--> spu_read_in_mbox(tid=%d)\n", tid);
  unsigned int cmd;

  pthread_mutex_lock(&msg_mutex[tid]);
  if (mbox_in[tid] == -1) {
    pthread_cond_wait(&msg_cond[tid], &msg_mutex[tid]);
  }
  cmd = mbox_in[tid];
  assert(cmd != -1);
  mbox_in[tid] = -1;
  pthread_mutex_unlock(&msg_mutex[tid]);
  //  printf("<-- spu_read_in_mbox() tid %d cmd = %#x\n", tid, cmd);
  return cmd;
}

int
spe_in_mbox_write(int tid, unsigned int *cmd, int l, int flags)
{
  //  printf("--> spe_in_mbox_write(%d, %d)\n", tid, *cmd);
  assert(l == 1);
  assert(flags == SPE_MBOX_ANY_NONBLOCKING);
  for (;;) {
    pthread_mutex_lock(&msg_mutex[tid]);
    if  (mbox_in[tid] == -1)
      break;
    pthread_mutex_unlock(&msg_mutex[tid]);
    sched_yield();
  }

  mbox_in[tid] = *cmd;
  pthread_cond_signal(&msg_cond[tid]);
  pthread_mutex_unlock(&msg_mutex[tid]);
  //  printf("<-- spe_in_mbox_write(%d, %d)\n", tid, *cmd);

  return l;
}

int
spe_out_mbox_status(int tid)
{
  int retval;

  pthread_mutex_lock(&mbox_out_mutex[tid]);
  if (mbox_out[tid] == -1) {
    retval = 0;
  } else {
    retval = 1;
  }
  pthread_mutex_unlock(&mbox_out_mutex[tid]);

  return retval;
}

int
spe_out_mbox_read(int tid, unsigned int *buf, int cnt)
{
  //  printf("--> spe_out_mbox_read(%d)\n", tid);
  for (;;) {
    pthread_mutex_lock(&mbox_out_mutex[tid]);
    if (mbox_out[tid] != -1)
      break;
    pthread_mutex_unlock(&mbox_out_mutex[tid]);
    sched_yield();
  }

  buf[0] = mbox_out[tid];
  mbox_out[tid] = -1;
  pthread_mutex_unlock(&mbox_out_mutex[tid]);

//  printf("<-- spe_out_mbox_read(%d, %d)\n", tid, buf[0]);
  return 1;
}

void
spu_write_out_mbox(unsigned int cmd)
{
  int tid = (int) (unsigned long) pthread_getspecific(key);
  //  printf("--> spu_write_out_mbox(tid=%d, %d)\n", tid, cmd);
  for (;;) {
    pthread_mutex_lock(&mbox_out_mutex[tid]);
    if (mbox_out[tid] == -1)
      break;

    pthread_mutex_unlock(&mbox_out_mutex[tid]);
    sched_yield();
  }
  mbox_out[tid] = cmd;
  pthread_mutex_unlock(&mbox_out_mutex[tid]);
  //  printf("<-- spu_write_out_mbox(tid=%d, %d)\n", tid, cmd);
}

spe_program_handle_t spe_progs[MAX_SPES];

int
spe_program_load(int speid, spe_program_handle_t *handle)
{
  spe_progs[speid] = *handle;
  
  return 0;
}

int
spe_context_run(int speid, unsigned int *entry, unsigned int runflags, void *argp,
		void *p1, void *p2)
{
  (spe_progs[speid])(speid, (unsigned long) argp, (unsigned long) p1);

  return 0;
}

int
spe_context_create(int i, void *p)
{
  static int __speid;
  int speid = __speid++;

  pthread_key_create(&key, NULL);
  pthread_setspecific(key, (void *)(unsigned long) speid);
  assert(speid < MAX_SPES);
  pthread_mutex_init(&msg_mutex[speid], NULL);
  pthread_cond_init (&msg_cond[speid], NULL);
  mbox_in[speid] = -1;
  mbox_out[speid] = -1;
  
  return speid;
}
