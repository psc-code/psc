
#include <mrc_bits.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>

static const int TERMS = 10; // no of terms to use in the Taylor series
 
static complex float
mrc_cerf(complex float z)
{
  const float cz = 2. / sqrt(M_PI);
  complex float res = 0., num, den;
  float fact = 1.;
  for(int i = 0; i < TERMS; i++) {
    num  = cz * pow(-1., i) * cpow(z, 2*i+1);
    den  = (2*i + 1) * fact;
    res  += num / den;
    fact *= i+1;
  }
 
  return res;
}

float
mrc_erfi(float x)
{
  return creal(-I * mrc_cerf(I * x));
}

void
mrc_cerf_test(void)
{
  printf( "Real error function\n\n");
  for ( float i = 0; i <= 1; i += 0.01 ){
    float gslerror = erf(i);
    float taylor   = mrc_cerf(i);
    printf("erf(%f): gsl = %f, taylor = %f, mag error = %f\n", i, gslerror, taylor,
	   fabs(gslerror-taylor));
  }
  
  complex float t, arg;
  printf( "\n\nImaginary error function\n\n");
  for (float i = 0; i <= 1; i += 0.01 ){
    arg = I * i;
    t   = mrc_cerf(arg);
    printf("erf(%f + i %f) = %f + i %f\n", creal(arg), cimag(arg), creal(t), cimag(t));
  }
}

