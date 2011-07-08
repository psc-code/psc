#include <stdarg.h>

static inline void clearfabi()
{
  FILE* f = fopen("fabi_out.txt", "wt");
  fclose(f);
}

static inline void printfabi(char *msg, ...)
{
  FILE* f=fopen("fabi_out.txt", "a+");
  va_list argList;
  va_start(argList, msg);

  vfprintf(f, msg, argList);

  va_end(argList);
  fclose(f);
  
  //vprintf(msg, argList);
}

#define HERE() {printfabi("[%s]\n", __func__);}

static inline void pause()
{
  while(getchar() != '\n');
}