
#ifndef FIELDS3D_HXX
#define FIELDS3D_HXX

template<typename R>
struct fields3d {
  using real_t = R;

  real_t *data;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  int first_comp; // first component
};


#endif

