#ifndef PSC_ALF_H
#define PSC_ALF_H
#include <stdio.h>
#include "psc_cbe.h"

#ifdef PSC_SPU_H
#include <alf_accel.h>
#else // PPU code
#include <alf.h>
void psc_alf_init_push_task_context(push_task_context_t * tcs, alf_task_desc_handle_t * task_desc, fields_cbe_t *pf);
///////
/// Access point for the shared ALF environment.
/// 
/// We don't want to pay the cost of starting up/stopping
/// the ALF environment for every task. Thus it's best to 
/// just set up a shared environment once, associate tasks
/// as needed, and destroy it at the end. It is initialized to
/// ALF null for the same reasons one initializes pointers 
/// to NULL.
extern alf_handle_t alf_env_shared;

#endif

#if CBE_DOUBLE
#define PSC_ALF_DATATYPE ALF_DATA_DOUBLE
#else
#define PSC_ALF_DATATYPE ALF_DATA_FLOAT
#endif 


/// Alf error checking
#define ACE { \
  if(ierr < 0){				\
  const char* err_str = alf_strerror(ierr);	\
  fprintf(stderr,"HERE: in %s() at %s:%d\n", __FUNCTION__, __FILE__, __LINE__); \
  fprintf(stderr,"ALF function returned an error: ");			\
  fprintf(stderr, err_str);						\
  fprintf(stderr, "\n");						\
  exit(2);								\
  }									\
  }

#endif
