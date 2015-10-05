
#ifndef MRC_BITS_H
#define MRC_BITS_H

#ifndef sqr
#define sqr(a) ((a) * (a))
#endif

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

float mrc_erfi(float x);

#endif
