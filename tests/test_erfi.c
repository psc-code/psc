
#include <math.h>
#include <mrc_bits.h>
#include <stdio.h>

int main(void)
{
  float x0 = .92, alpha = 1.85;
  int N = 1000;
  for (int i = 0; i <= N; i++) {
    float x = -1.5 + 3./N * i;
    float B;
    if (x < -x0) {
      B = -1.;
    } else if (x > x0) {
      B = 1.;
    } else {
      B = alpha * exp(-sqr(x)) * sqrt(M_PI) / 2. * mrc_erfi(x);
    }
    printf("%g %g\n", x, B);
  }

  return 0;
}
