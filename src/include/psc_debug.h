
#ifndef PSC_DEBUG_H
#define PSC_DEBUG_H

static inline void
float_set_nan(float *x)
{
  *(int *) x = 0x7fbfffff;
}

static inline void
double_set_nan(double *x)
{
  *(long long *) x = 0x7ff8000000000000;
}

#endif
