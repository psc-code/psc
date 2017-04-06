
#ifndef PDE_FORTRAN_H
#define PDE_FORTRAN_H

#ifndef FORTRAN_TYPEDEFS
#define FORTRAN_TYPEDEFS

typedef float real;
typedef int integer;
typedef int logical;

#endif

#define F(p_f, m) (&F3S(p_f, m, -s_sw[0], -s_sw[1], -s_sw[2]))

#endif

