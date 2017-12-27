
#ifndef PSC_HYDRO_ARRAY_H
#define PSC_HYDRO_ARRAY_H

// Generic looping
#define XYZ_LOOP(xl,xh,yl,yh,zl,zh) \
  for( z=zl; z<=zh; z++ )	    \
    for( y=yl; y<=yh; y++ )	    \
      for( x=xl; x<=xh; x++ )
	      
// x_NODE_LOOP => Loop over all non-ghost nodes at plane x
#define x_NODE_LOOP(x) XYZ_LOOP(x,x,1,ny+1,1,nz+1)
#define y_NODE_LOOP(y) XYZ_LOOP(1,nx+1,y,y,1,nz+1)
#define z_NODE_LOOP(z) XYZ_LOOP(1,nx+1,1,ny+1,z,z)

// ======================================================================
// PscHydroArray

template<class HydroArrayBase>
struct PscHydroArray : HydroArrayBase
{
  typedef HydroArrayBase Base;
  using typename Base::Grid;
  using typename Base::Element;

  using Base::data;
  using Base::grid;

  static PscHydroArray* create(Grid *grid)
  {
    return static_cast<PscHydroArray*>(Base::create(grid));
  }

  // ----------------------------------------------------------------------
  // clear
  
  void clear()
  {
    Element *h = data();
    memset(h, 0, grid()->nv * sizeof(*h));
  }

  // ----------------------------------------------------------------------
  // synchronize
  
  template<class G, class F3D>
  struct CommHydro : Comm<G, F3D>
  {
    typedef Comm<G, F3D> Base;
    typedef G Grid;
    
    using Base::begin;
    using Base::end;

    using Base::nx_;
    using Base::g_;
    using Base::buf_size_;

    CommHydro(Grid *g) : Base(g)
    {
      for (int X = 0; X < 3; X++) {
	int Y = (X + 1) % 3, Z = (X + 2) % 3;
	buf_size_[X] = 14 * (nx_[Y] + 1) * (nx_[Z] + 1);
      }
    }

    void begin_send(int X, int side, float* p, F3D& F)
    {
      int face = side ? nx_[X] + 1 : 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) {
	  Element* h = &F(x,y,z);
	  *p++ = h->jx;
	  *p++ = h->jy;
	  *p++ = h->jz;
	  *p++ = h->rho;
	  *p++ = h->px;
	  *p++ = h->py;
	  *p++ = h->pz;
	  *p++ = h->ke;
	  *p++ = h->txx;
	  *p++ = h->tyy;
	  *p++ = h->tzz;
	  *p++ = h->tyz;
	  *p++ = h->tzx;
	  *p++ = h->txy;
	});
    }
    
    void end_recv(int X, int side, float* p, F3D& F)
    {
      int face = side ? 1 : nx_[X] + 1;
      foreach_node(g_, X, face, [&](int x, int y, int z) {
	  Element* h = &F(x,y,z);
	  h->jx  += *p++;
	  h->jy  += *p++;
	  h->jz  += *p++;
	  h->rho += *p++;
	  h->px  += *p++;
	  h->py  += *p++;
	  h->pz  += *p++;
	  h->ke  += *p++;
	  h->txx += *p++;
	  h->tyy += *p++;
	  h->tzz += *p++;
	  h->tyz += *p++;
	  h->tzx += *p++;
	  h->txy += *p++;
	});
    }
  };

  void synchronize()
  {
    int face, bc, x, y, z;
    Element* h;

    const Grid* g = grid();
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    Field3D<HydroArrayBase> H(*this);

    // Note: synchronize_hydro assumes that hydro has not been adjusted
    // at the local domain boundary. Because hydro fields are purely
    // diagnostic, correct the hydro along local boundaries to account
    // for accumulations over partial cell volumes

# define ADJUST_HYDRO(i,j,k,X,Y,Z)              \
    do {					\
      bc = g->bc[BOUNDARY(i,j,k)];		\
      if( bc<0 || bc>=psc_world_size ) {	\
	face = (i+j+k)<0 ? 1 : n##X+1;		\
	X##_NODE_LOOP(face) {			\
	  h = &H(x,y,z);			\
	  h->jx  *= 2;				\
	  h->jy  *= 2;				\
	  h->jz  *= 2;				\
	  h->rho *= 2;				\
	  h->px  *= 2;				\
	  h->py  *= 2;				\
	  h->pz  *= 2;				\
	  h->ke  *= 2;				\
	  h->txx *= 2;				\
	  h->tyy *= 2;				\
	  h->tzz *= 2;				\
	  h->tyz *= 2;				\
	  h->tzx *= 2;				\
	  h->txy *= 2;				\
	}					\
      }						\
    } while(0)
  
    ADJUST_HYDRO(-1, 0, 0,x,y,z);
    ADJUST_HYDRO( 0,-1, 0,y,z,x);
    ADJUST_HYDRO( 0, 0,-1,z,x,y);
    ADJUST_HYDRO( 1, 0, 0,x,y,z);
    ADJUST_HYDRO( 0, 1, 0,y,z,x);
    ADJUST_HYDRO( 0, 0, 1,z,x,y);

# undef ADJUST_HYDRO

    CommHydro<Grid, Field3D<HydroArrayBase>> comm(this->grid());

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, H);
      comm.end(dir, H);
    }
    
  }
};

#endif

