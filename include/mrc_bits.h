
#ifndef MRC_BITS_H
#define MRC_BITS_H

#define sqr(a) ((a) * (a))

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

float mrc_erfi(float x);

#endif
