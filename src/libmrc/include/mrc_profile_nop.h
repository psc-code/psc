
#ifndef PROFILE_NOP_H
#define PROFILE_NOP_H

// provides to profile.h API, but gives only nop implementations.
// to turn off profiling in certain files, include this file instead of
// the regular profile.h

static inline int
prof_register(const char *name, float simd, int flops, int bytes)
{
  return 1;
}

static inline void
prof_start(int pr)
{
}

static inline void
prof_stop(int pr)
{
}


#endif
