
#ifndef GRID_H
#define GRID_H

#define private public
#include <vpic.h>

// ======================================================================
// Grid

struct Grid {
  Grid(grid_t *grid);
  
  void setup(double dx[3], double dt, double cvac, double eps0);
  void partition_periodic_box(double xl[3], double xh[3], int gdims[3], int np[3]);
  void set_fbc(int boundary, int fbc);
  void set_pbc(int boundary, int pbc);
  void mp_size_recv_buffer(int tag, int size);
  void mp_size_send_buffer(int tag, int size);

  grid_t* getGrid_t() { return g_; }

private:
  grid_t *g_;
};

inline Grid::Grid(grid_t *g)
  : g_(g)
{
}

inline void Grid::setup(double dx[3], double dt, double cvac, double eps0)
{
  g_->dx = dx[0];
  g_->dy = dx[1];
  g_->dz = dx[2];
  g_->dt = dt;
  g_->cvac = cvac;
  g_->eps0 = eps0;
}

inline void Grid::partition_periodic_box(double xl[3], double xh[3],
					 int gdims[3], int np[3])
{
  ::partition_periodic_box(g_, xl[0], xl[1], xl[2], xh[0], xh[1], xh[2],
			   gdims[0], gdims[1], gdims[2], np[0], np[1], np[2]);
}

inline void Grid::set_fbc(int boundary, int fbc)
{
  ::set_fbc(g_, boundary, fbc);
}

inline void Grid::set_pbc(int boundary, int pbc)
{
  ::set_pbc(g_, boundary, pbc);
}

inline void Grid::mp_size_recv_buffer(int tag, int size)
{
  ::mp_size_recv_buffer(g_->mp, tag, size);
}

inline void Grid::mp_size_send_buffer(int tag, int size)
{
  ::mp_size_send_buffer(g_->mp, tag, size);
}

#endif
