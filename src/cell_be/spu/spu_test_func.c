#include <stdio.h> 
#include "spu_test_func.h"

int spu_psc_hello(void) 
{ 
  printf("%s from [%#llx]\n",global_ctx.hello,global_ctx.spe_id);
  return 0;
}

int spu_psc_goodbye(void) 
{
  printf("%s from [%#llx]\n",global_ctx.bye,global_ctx.spe_id);
  return 0;
}


