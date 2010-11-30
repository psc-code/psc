#include <stdio.h> 
#include "spu_test_func.h"

int spu_psc_hello(void) 
{ 
  printf("%s from [%#llx]\n",psc_env.hello,psc_env.spe_id);
  return 0;
}

int spu_psc_goodbye(void) 
{
  printf("%s from [%#llx]\n",psc_env.bye,psc_env.spe_id);
  return 0;
}


