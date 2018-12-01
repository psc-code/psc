
#pragma once

#include "GridLoop.h"
#include "Field3D.h"

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
// PscHydroyOps

template<typename _Mparticles, typename _MfieldsHydro, typename _MfieldsInterpolator>
struct PscHydroOps
{
  using Mparticles = _Mparticles;
  using MfieldsHydro = _MfieldsHydro;
  using MfieldsInterpolator = _MfieldsInterpolator;
  using Grid = typename MfieldsHydro::Grid;
  using Element = typename MfieldsHydro::Element;
  
  // ----------------------------------------------------------------------
  // clear
  
  static void clear(MfieldsHydro& hydro)
  {
    auto h = hydro.data();
    memset(h, 0, hydro.vgrid()->nv * sizeof(*h) * MfieldsHydro::N_COMP);
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

    CommHydro(Grid& g) : Base(g)
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

  static void synchronize(MfieldsHydro& hydro)
  {
    using F3D = Field3D<typename MfieldsHydro::Patch>;
    int face, bc, x, y, z;
    Element* h;

    auto g = hydro.vgrid();
    const int nx = g->nx, ny = g->ny, nz = g->nz;
    auto fields_hydro = hydro[0];
    F3D H{hydro.getPatch(0)};

    // Note: synchronize_hydro assumes that hydro has not been adjusted
    // at the local domain boundary. Because hydro fields are purely
    // diagnostic, correct the hydro along local boundaries to account
    // for accumulations over partial cell volumes

# define ADJUST_HYDRO(i,j,k,X,Y,Z)              \
    do {					\
      int bc = g->bc[BOUNDARY(i,j,k)];		\
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

    CommHydro<Grid, F3D> comm{*hydro.vgrid()};

    for (int dir = 0; dir < 3; dir++) {
      comm.begin(dir, H);
      comm.end(dir, H);
    }
  }

  // ----------------------------------------------------------------------
  // accumulate_hydro_p
  
  // accumulate_hydro_p adds the hydrodynamic fields associated with the
  // supplied particle_list to the hydro array.  Trilinear interpolation
  // is used.  hydro is known at the nodes at the same time as particle
  // positions. No effort is made to fix up edges of the computational
  // domain.  All particles on the list must be inbounds.  Note, the
  // hydro jx,jy,jz are for diagnostic purposes only; they are not
  // accumulated with a charge conserving algorithm.

  static void accumulate_hydro_p(MfieldsHydro& hydro, const typename Mparticles::const_iterator sp_iter,
				 /*const*/ MfieldsInterpolator& interpolator)
  {
    auto& sp = *sp_iter;
    auto& ha = hydro.getPatch(0);
    auto& IP = interpolator.getPatch(0);
    float c, qsp, mspc, qdt_2mc, qdt_4mc2, r8V;
    int np, stride_10, stride_21, stride_43;

    float dx, dy, dz, ux, uy, uz, w, vx, vy, vz, ke_mc;
    float w0, w1, w2, w3, w4, w5, w6, w7, t;
    int i, n;

    const auto& g = sp.vgrid();
    const auto* p = sp.p;

    c        = g.cvac;
    qsp      = sp.q;
    mspc     = sp.m*c;
    qdt_2mc  = (qsp*g.dt)/(2*mspc);
    qdt_4mc2 = qdt_2mc / (2*c);
    r8V      = g.r8V;

    np        = sp.np;
    stride_10 = (VOXEL(1,0,0, g.nx,g.ny,g.nz) -
		 VOXEL(0,0,0, g.nx,g.ny,g.nz));
    stride_21 = (VOXEL(0,1,0, g.nx,g.ny,g.nz) -
		 VOXEL(1,0,0, g.nx,g.ny,g.nz));
    stride_43 = (VOXEL(0,0,1, g.nx,g.ny,g.nz) -
		 VOXEL(1,1,0, g.nx,g.ny,g.nz));

    for( n=0; n<np; n++ ) {

      // Load the particle
      dx = p[n].dx;
      dy = p[n].dy;
      dz = p[n].dz;
      i  = p[n].i;
      ux = p[n].ux;
      uy = p[n].uy;
      uz = p[n].uz;
      w  = p[n].w;
    
      // Half advance E
      ux += qdt_2mc*((IP[i].ex+dy*IP[i].dexdy) + dz*(IP[i].dexdz+dy*IP[i].d2exdydz));
      uy += qdt_2mc*((IP[i].ey+dz*IP[i].deydz) + dx*(IP[i].deydx+dz*IP[i].d2eydzdx));
      uz += qdt_2mc*((IP[i].ez+dx*IP[i].dezdx) + dy*(IP[i].dezdy+dx*IP[i].d2ezdxdy));

      // Boris rotation - Interpolate B field
      w5 = IP[i].cbx + dx*IP[i].dcbxdx;
      w6 = IP[i].cby + dy*IP[i].dcbydy;
      w7 = IP[i].cbz + dz*IP[i].dcbzdz;

      // Boris rotation - curl scalars (0.5 in v0 for half rotate) and
      // kinetic energy computation. Note: gamma-1 = |u|^2 / (gamma+1)
      // is the numerically accurate way to compute gamma-1
      ke_mc = ux*ux + uy*uy + uz*uz; // ke_mc = |u|^2 (invariant)
      vz = sqrt(1+ke_mc);            // vz = gamma    (invariant)
      ke_mc *= c/(vz+1);             // ke_mc = c|u|^2/(gamma+1) = c*(gamma-1)
      vz = c/vz;                     // vz = c/gamma
      w0 = qdt_4mc2*vz;
      w1 = w5*w5 + w6*w6 + w7*w7;    // |cB|^2
      w2 = w0*w0*w1;
      w3 = w0*(1+(1./3.)*w2*(1+0.4*w2));
      w4 = w3/(1 + w1*w3*w3); w4 += w4;

      // Boris rotation - uprime
      w0 = ux + w3*( uy*w7 - uz*w6 );
      w1 = uy + w3*( uz*w5 - ux*w7 );
      w2 = uz + w3*( ux*w6 - uy*w5 );

      // Boris rotation - u
      ux += w4*( w1*w7 - w2*w6 );
      uy += w4*( w2*w5 - w0*w7 );
      uz += w4*( w0*w6 - w1*w5 );

      // Compute physical velocities
      vx  = ux*vz;
      vy  = uy*vz;
      vz *= uz;

      // Compute the trilinear coefficients
      w0  = r8V*w;    // w0 = (1/8)(w/V)
      dx *= w0;       // dx = (1/8)(w/V) x
      w1  = w0+dx;    // w1 = (1/8)(w/V) + (1/8)(w/V)x = (1/8)(w/V)(1+x)
      w0 -= dx;       // w0 = (1/8)(w/V) - (1/8)(w/V)x = (1/8)(w/V)(1-x)
      w3  = 1+dy;     // w3 = 1+y
      w2  = w0*w3;    // w2 = (1/8)(w/V)(1-x)(1+y)
      w3 *= w1;       // w3 = (1/8)(w/V)(1+x)(1+y)
      dy  = 1-dy;     // dy = 1-y
      w0 *= dy;       // w0 = (1/8)(w/V)(1-x)(1-y)
      w1 *= dy;       // w1 = (1/8)(w/V)(1+x)(1-y)
      w7  = 1+dz;     // w7 = 1+z
      w4  = w0*w7;    // w4 = (1/8)(w/V)(1-x)(1-y)(1+z) = (w/V) trilin_0 *Done
      w5  = w1*w7;    // w5 = (1/8)(w/V)(1+x)(1-y)(1+z) = (w/V) trilin_1 *Done
      w6  = w2*w7;    // w6 = (1/8)(w/V)(1-x)(1+y)(1+z) = (w/V) trilin_2 *Done
      w7 *= w3;       // w7 = (1/8)(w/V)(1+x)(1+y)(1+z) = (w/V) trilin_3 *Done
      dz  = 1-dz;     // dz = 1-z
      w0 *= dz;       // w0 = (1/8)(w/V)(1-x)(1-y)(1-z) = (w/V) trilin_4 *Done
      w1 *= dz;       // w1 = (1/8)(w/V)(1+x)(1-y)(1-z) = (w/V) trilin_5 *Done
      w2 *= dz;       // w2 = (1/8)(w/V)(1-x)(1+y)(1-z) = (w/V) trilin_6 *Done
      w3 *= dz;       // w3 = (1/8)(w/V)(1+x)(1+y)(1-z) = (w/V) trilin_7 *Done

      // Accumulate the hydro fields
#   define ACCUM_HYDRO( wn)						\
      t  = qsp*wn;        /* t  = (qsp w/V) trilin_n */			\
      ha[i].jx  += t*vx;						\
      ha[i].jy  += t*vy;						\
      ha[i].jz  += t*vz;						\
      ha[i].rho += t;							\
      t  = mspc*wn;       /* t = (msp c w/V) trilin_n */		\
      dx = t*ux;          /* dx = (px w/V) trilin_n */			\
      dy = t*uy;							\
      dz = t*uz;							\
      ha[i].px  += dx;							\
      ha[i].py  += dy;							\
      ha[i].pz  += dz;							\
      ha[i].ke  += t*ke_mc;						\
      ha[i].txx += dx*vx;						\
      ha[i].tyy += dy*vy;						\
      ha[i].tzz += dz*vz;						\
      ha[i].tyz += dy*vz;						\
      ha[i].tzx += dz*vx;						\
      ha[i].txy += dx*vy

      /**/            ACCUM_HYDRO(w0); // Cell i,j,k
      i += stride_10; ACCUM_HYDRO(w1); // Cell i+1,j,k
      i += stride_21; ACCUM_HYDRO(w2); // Cell i,j+1,k
      i += stride_10; ACCUM_HYDRO(w3); // Cell i+1,j+1,k
      i += stride_43; ACCUM_HYDRO(w4); // Cell i,j,k+1
      i += stride_10; ACCUM_HYDRO(w5); // Cell i+1,j,k+1
      i += stride_21; ACCUM_HYDRO(w6); // Cell i,j+1,k+1
      i += stride_10; ACCUM_HYDRO(w7); // Cell i+1,j+1,k+1

#   undef ACCUM_HYDRO
    }
  }
};

#undef XYZ_LOOP
#undef x_NODE_LOOP
#undef y_NODE_LOOP
#undef z_NODE_LOOP

