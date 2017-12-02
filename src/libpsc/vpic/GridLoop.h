
#ifndef GRID_LOOP_H
#define GRID_LOOP_H

// ======================================================================
// foreach

template<class F>
static void foreach(F f, int ib, int ie, int jb, int je, int kb, int ke)
{
  for (int k = kb; k <= ke; k++) {
    for (int j = jb; j <= je; j++) {
      for (int i = ib; i <= ie; i++) {
	f(i,j,k);
      }
    }
  }
};

template<class Grid, class F>
static void foreach_edge(const Grid *g, int Y, int Z, int face, F f)
{
  if        (Y == 1 && Z == 2) { foreach(f, face, face, 1, g->ny  , 1, g->nz+1);
  } else if (Y == 2 && Z == 1) { foreach(f, face, face, 1, g->ny+1, 1, g->nz  );
  } else if (Y == 2 && Z == 0) { foreach(f, 1, g->nx+1, face, face, 1, g->nz  );
  } else if (Y == 0 && Z == 2) { foreach(f, 1, g->nx  , face, face, 1, g->nz+1);
  } else if (Y == 0 && Z == 1) { foreach(f, 1, g->nx  , 1, g->ny+1, face, face);
  } else if (Y == 1 && Z == 0) { foreach(f, 1, g->nx+1, 1, g->ny  , face, face);
  } else {
    assert(0);
  }
}
  
template<class Grid, class F>
static void foreach_node(const Grid *g, int X, int face, F f)
{
  if        (X == 0) { foreach(f, face, face, 1, g->ny+1, 1, g->nz+1);
  } else if (X == 1) { foreach(f, 1, g->nx+1, face, face, 1, g->nz+1);
  } else if (X == 2) { foreach(f, 1, g->nx+1, 1, g->ny+1, face, face);
  } else {
    assert(0);
  }
}

template<class Grid, class F>
static void foreach_face(const Grid *g, int X, int face, F f)
{
  if        (X == 0) { foreach(f, face, face, 1, g->ny  , 1, g->nz  );
  } else if (X == 1) { foreach(f, 1, g->nx  , face, face, 1, g->nz  );
  } else if (X == 2) { foreach(f, 1, g->nx  , 1, g->ny  , face, face);
  } else {
    assert(0);
  }
}

template<class Grid, class F>
static void foreach_nc_interior(F f, Grid* g)
{
  const int nx = g->nx, ny = g->ny, nz = g->nz;
  foreach(f, 2, nx, 2, ny, 2, nz);
}

template<class Grid, class F>
static void foreach_nc_boundary(F f, Grid* g)
{
  const int nx = g->nx, ny = g->ny, nz = g->nz;

  // z faces, x edges, y edges and all corners
  foreach(f, 1   , nx+1, 1   , ny+1, 1   , 1   );
  foreach(f, 1   , nx+1, 1   , ny+1, nz+1, nz+1);

  // y faces, z edges
  foreach(f, 1   , nx+1, 1   , 1   ,    2, nz  );
  foreach(f, 1   , nx+1, ny+1, ny+1,    2, nz  );

  // x faces
  foreach(f, 1   , 1   , 2   , ny  ,    2, nz  );
  foreach(f, nx+1, nx+1, 2   , ny  ,    2, nz  );
}

template<class Grid, class F>
static void foreach_ec_interior(F f, Grid* g)
{
  const int nx = g->nx, ny = g->ny, nz = g->nz;

  for (int k = 2; k <= nz; k++) {
    for (int j = 2; j <= ny; j++) {
      for (int i = 2; i <= nx; i++) {
	f.x(i,j,k);
	f.y(i,j,k);
	f.z(i,j,k);
      }
    }
  }

  // Do leftover interior ex
  for (int k = 2; k <= nz; k++) {
    for (int j = 2; j <= ny; j++) {
      f.x(1,j,k);
    }
  }

  // Do leftover interior ey
  for (int k = 2; k <= nz; k++) {
    for (int i = 2; i <= nx; i++) {
      f.y(i,1,k);
    }
  }

  // Do leftover interior ez
  for (int j = 2; j <= ny; j++) {
    for (int i = 2; i <= nx; i++) {
      f.z(i,j,1);
    }
  }
}

template<class Grid, class F>
static void foreach_ec_boundary(F f, Grid* g)
{
  const int nx = g->nx, ny = g->ny, nz = g->nz;

  for (int j = 1; j <= ny+1; j++) {
    for (int i = 1; i <= nx; i++) {
      f.x(i,j,1);
    }
  }
  for (int j = 1; j <= ny+1; j++) {
    for (int i = 1; i <= nx; i++) {
      f.x(i,j,nz+1);
    }
  }
  for (int k = 2; k <= nz; k++) {
    for (int i = 1; i <= nx; i++) {
      f.x(i,1,k);
    }
  }
  for (int k = 2; k <= nz; k++) {
    for (int i = 1; i <= nx; i++) {
      f.x(i,ny+1,k);
    }
  }

  // Do exterior ey
  for (int k = 1; k <= nz+1; k++) {
    for (int j = 1; j <= ny; j++) {
      f.y(1,j,k);
    }
  }
  for (int k = 1; k <= nz+1; k++) {
    for (int j = 1; j <= ny; j++) {
      f.y(nx+1,j,k);
    }
  }
  for (int j = 1; j <= ny; j++) {
    for (int i = 2; i <= nx; i++) {
      f.y(i,j,1);
    }
  }
  for (int j = 1; j <= ny; j++) {
    for (int i = 2; i <= nx; i++) {
      f.y(i,j,nz+1);
    }
  }

  // Do exterior ez
  for (int k = 1; k <= nz; k++) {
    for (int i = 1; i <= nx+1; i++) {
      f.z(i,1,k);
    }
  }
  for (int k = 1; k <= nz; k++) {
    for (int i = 1; i <= nx+1; i++) {
      f.z(i,ny+1,k);
    }
  }
  for (int k = 1; k <= nz; k++) {
    for (int j = 2; j <= ny; j++) {
      f.z(1,j,k);
    }
  }
  for (int k = 1; k <= nz; k++) {
    for (int j = 2; j <= ny; j++) {
      f.z(nx+1,j,k);
    }
  }
}
  

// ======================================================================
// Comm

template<class G, class F3D>
struct Comm
{
  typedef G Grid;
  
  Comm(Grid* g) : g_(g)
  {
    nx_[0] = g_->nx;
    nx_[1] = g_->ny;
    nx_[2] = g_->nz;
  }

  const int to_ijk[3][2][3] = {
    { { -1, 0, 0 },
      {  1, 0, 0 }, },
    { {  0,-1, 0 },
      {  0, 1, 0 }, },
    { {  0, 0,-1 },
      {  0, 0, 1 }, },
  };
  // wrap what probably should be in Grid::

  float* get_send_buf(int dir, int side, int sz) const
  {
    const int* ijk = to_ijk[dir][side];
    return static_cast<float*>(g_->size_send_port(ijk[0], ijk[1], ijk[2], sz * sizeof(float)));
  }

  void begin_send_port(int dir, int side, int sz) const
  {
    const int* ijk = to_ijk[dir][side];
    g_->begin_send_port(ijk[0], ijk[1], ijk[2], sz * sizeof(float));
  }

  void begin_recv_port(int dir, int side, int sz) const
  {
    const int* ijk = to_ijk[dir][side];
    g_->begin_recv_port(ijk[0], ijk[1], ijk[2], sz * sizeof(float));
  }

  float* end_recv_port(int dir, int side) const
  {
    const int* ijk = to_ijk[dir][side];
    return static_cast<float*>(g_->end_recv_port(ijk[0], ijk[1], ijk[2]));
  }

  void end_send_port(int dir, int side) const
  {
    const int* ijk = to_ijk[dir][side];
    g_->end_send_port(ijk[0], ijk[1], ijk[2]);
  }

  // ----------------------------------------------------------------------

  virtual void begin_send(int dir, int side, float* p, F3D& F) = 0;
  virtual void end_recv(int dir, int side, float* p, F3D& F) = 0;

  void begin_recv(int dir, int side)
  {
    begin_recv_port(dir, side, buf_size_[dir]);
  }

  void begin_send(int dir, int side, F3D& F)
  {
    float *p = get_send_buf(dir, side, buf_size_[dir]);
    if (p) {
      begin_send(dir, side, p, F);
      begin_send_port(dir, side, buf_size_[dir]);
    }
  }

  void end_recv(int dir, int side, F3D& F)
  {
    float* p = end_recv_port(dir, side);
    if (p) {
      end_recv(dir, side, p, F);
    }
  }
  
  void end_send(int dir, int side)
  {
    end_send_port(dir, side);
  }
  
  void begin(int dir, F3D& F)
  {
    for (int side = 0; side <= 1; side++) {
      begin_recv(dir, side);
    }
    for (int side = 0; side <= 1; side++) {
      begin_send(dir, side, F);
    }
  }

  void end(int dir, F3D& F)
  {
    for (int side = 0; side <= 1; side++) {
      end_recv(dir, side, F);
    }

    for (int side = 0; side <= 1; side++) {
      end_send(dir, side);
    }
  }
  
  void begin(F3D& F)
  {
    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	begin_recv(dir, side);
      }
    }

    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	begin_send(dir, side, F);
      }
    }
  }
  
  void end(F3D& F)
  {
    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	end_recv(dir, side, F);
      }
    }

    for (int side = 0; side <= 1; side++) {
      for (int dir = 0; dir < 3; dir++) {
	end_send(dir, side);
      }
    }
  }
  
protected:
  int nx_[3];
  int buf_size_[3];
  Grid *g_;
};

#endif

