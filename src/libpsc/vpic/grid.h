
#ifndef GRID_H
#define GRID_H

#include <vpic.h>

// ======================================================================
// Grid

struct Grid : grid_t {
  void setup(double dx[3], double dt, double cvac, double eps0);
  void partition_periodic_box(double xl[3], double xh[3], int gdims[3], int np[3]);
  void set_fbc(int boundary, int fbc);
  void set_pbc(int boundary, int pbc);
  void mp_size_recv_buffer(int tag, int size);
  void mp_size_send_buffer(int tag, int size);

  void* size_send_port(int i, int j, int k, int sz)
  {
    return ::size_send_port(i, j, k, sz, this);
  }
  
  void begin_send_port(int i, int j, int k, int sz)
  {
    ::begin_send_port(i, j, k, sz, this);
  }

  void end_send_port(int i, int j, int k)
  {
    ::end_send_port(i, j, k, this);
  }

  void begin_recv_port(int i, int j, int k, int sz)
  {
    ::begin_recv_port(i, j, k, sz, this);
  }

  void* end_recv_port(int i, int j, int k)
  {
    return ::end_recv_port(i, j, k, this);
  }
  
};

inline void Grid::setup(double dx_[3], double dt_, double cvac_, double eps0_)
{
  dx = dx_[0];
  dy = dx_[1];
  dz = dx_[2];
  dt = dt_;
  cvac = cvac_;
  eps0 = eps0_;
}

inline void Grid::partition_periodic_box(double xl[3], double xh[3],
					 int gdims[3], int np[3])
{
  ::partition_periodic_box(this, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			   gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
}

inline void Grid::set_fbc(int boundary, int fbc)
{
  ::set_fbc(this, boundary, fbc);
}

inline void Grid::set_pbc(int boundary, int pbc)
{
  ::set_pbc(this, boundary, pbc);
}

inline void Grid::mp_size_recv_buffer(int tag, int size)
{
  ::mp_size_recv_buffer(mp, tag, size);
}

inline void Grid::mp_size_send_buffer(int tag, int size)
{
  ::mp_size_send_buffer(mp, tag, size);
}

#endif
