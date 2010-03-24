//#include <iostream>
// kay
#include <cmath>

extern "C"
{
  void setpotential_(double *pot, double *x, double *y, double *z, double *t,
		     double *ld, double *wl)
  {
    double xPos = 5.0e-6/ *ld;
    double yPos = 5.0e-6/ *ld;
    double zPos1 = 0.0e-6/ *ld;
    double zPos2 = 10.0e-6/ *ld;
    double dist = 1.0e-7/ *ld;
    double omega = *wl;
//     if (fabs(xPos - *x) < dist) {
//       if (fabs(yPos - *y) < dist) {
 	if (fabs(zPos1 - *z) < dist) {
 	  *pot = 0.0;
	}
 	if (fabs(zPos2 - *z) < dist) {
 	  *pot = 1.0*sin(omega * *t / *wl);
	}
	//      }
	//    }
  }
}
