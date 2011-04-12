#include <stdio.h> 
#include "psc_spu.h"
#include "spu_test_func.h"

int spu_psc_hello(void) 
{ 
  fprintf(stdout,"%s from [%#llx]\n",spu_ctx.hello,spu_ctx.spe_id);
  fflush(stdout);
  return 0;
}

int spu_psc_goodbye(void) 
{
  fprintf(stdout,"%s from [%#llx]\n",spu_ctx.bye,spu_ctx.spe_id);
  fflush(stdout);
  return 0;
}


