//#include <iostream>
// kay
#include <cmath>

extern "C"
{
  void fdtdsolver_(double *pot0, double *pot1, double *rho,
		   int *xMin, int *xMax,
		   int *yMin, int *yMax, int *zMin, int *zMax,
		   double *cnx, double *cny, double *cnz, double *dt)
  {
    int dimX = *xMax-*xMin+1;
    int dimY = *yMax-*yMin+1;
    int dimZ = *zMax-*zMin+1;
    double tmp;
    double ctx = *cnx * *cnx;
    double cty = *cny * *cny;
    double ctz = *cnz * *cnz;
    double ctxyz2 = 2 - 2*(ctx + cty + ctz);
    double ctt = *dt * *dt * 0.5;

//////// fdtdloop
    for(int l=1; l<dimZ-1; ++l) {
      for(int k=1; k<dimY-1; ++k) {
	for(int j=1; j<dimX-1; ++j) {
	  pot0[j+dimX*(k+dimY*l)] =
	    ctz*(pot1[j+dimX*(k+dimY* (l-1) )]+pot1[j+dimX*(k+dimY* (l+1) )])
	    + cty*(pot1[j+dimX*( k-1 +dimY*l)]+pot1[j+dimX*( k+1 +dimY*l)])
	    + ctx*(pot1[ j-1 +dimX*(k+dimY*l)]+pot1[ j+1 +dimX*(k+dimY*l)])
	    + ctxyz2*pot1[j+dimX*(k+dimY*l)]-pot0[j+dimX*(k+dimY*l)]
	    + ctt*rho[j+dimX*(k+dimY*l)];
	}
      }
    }
//////// swap arrays
    for(int l=1; l<dimZ-1; ++l) {
      for(int k=1; k<dimY-1; ++k) {
	for(int j=1; j<dimX-1; ++j) {
	  tmp = pot0[j+dimX*(k+dimY*l)];
	  pot0[j+dimX*(k+dimY*l)] = pot1[j+dimX*(k+dimY*l)];
	  pot1[j+dimX*(k+dimY*l)] = tmp;
	}
      }
    }
//     for(int l=1; l<dimZ-1; ++l) {
//       for(int k=1; k<dimY-1; ++k) {
// 	for(int j=1; j<dimX-1; ++j) {
// 	  pot1[j+dimX*(k+dimY*l)] =
// 	    ctx*(pot0[j+dimX*(k+dimY* (l-1) )]+pot0[j+dimX*(k+dimY* (l+1) )])
// 	    + cty*(pot0[j+dimX*( k-1 +dimY*l)]+pot0[j+dimX*( k+1 +dimY*l)])
// 	    + ctz*(pot0[ j-1 +dimX*(k+dimY*l)]+pot0[ j+1 +dimX*(k+dimY*l)])
// 	    + ctxyz2*pot0[j+dimX*(k+dimY*l)]-pot1[j+dimX*(k+dimY*l)];
	  
// 	}
//       }
//     }
  }
}
