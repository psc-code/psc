
#ifndef PSC_PARTICLES_BASE_H
#define PSC_PARTICLES_BASE_H

#include "VpicListBase.h"

// ======================================================================
// PscSpecies

// FIXME: Eventually particle_t (definitely) and the other formats
// (maybe) should be opaque and specific to a particular
// species_advance implementation

struct PscParticle
{
  float dx, dy, dz; // Particle position in cell coordinates (on [-1,1])
  int32_t i;        // Voxel containing the particle.  Note that
  /**/              // particled awaiting processing by boundary_p
  /**/              // have actually set this to 8*voxel + face where
  /**/              // face is the index of the face they interacted
  /**/              // with (on 0:5).  This limits the local number of
  /**/              // voxels to 2^28 but emitter handling already
  /**/              // has a stricter limit on this (2^26).
  float ux, uy, uz; // Particle normalized momentum
  float w;          // Particle weight (number of physical particles)
};

// WARNING: FUNCTIONS THAT USE A PARTICLE_MOVER ASSUME THAT EVERYBODY
// WHO USES THAT PARTICLE MOVER WILL HAVE ACCESS TO PARTICLE ARRAY

struct PscParticleMover
{
  float dispx, dispy, dispz; // Displacement of particle
  int32_t i;                 // Index of the particle to move
};

// NOTE: THE LAYOUT OF A PARTICLE_INJECTOR _MUST_ BE COMPATIBLE WITH
// THE CONCATENATION OF A PARTICLE_T AND A PARTICLE_MOVER!

struct PscParticleInjector
{
  float dx, dy, dz;          // Particle position in cell coords (on [-1,1])
  int32_t i;                 // Index of cell containing the particle
  float ux, uy, uz;          // Particle normalized momentum
  float w;                   // Particle weight (number of physical particles)
  float dispx, dispy, dispz; // Displacement of particle
  SpeciesId sp_id;           // Species of particle
};

template<class G>
struct PscSpecies
{
  typedef G Grid;
  
  // ----------------------------------------------------------------------
  // PscSpecies::ctor
  
  PscSpecies(const char * name, float q, float m,
	     int max_local_np, int max_local_nm,
	     int sort_interval, int sort_out_of_place,
	     Grid *grid)
  {
    memset(this, 0, sizeof(*this)); // FIXME
    int len = name ? strlen(name) : 0;
    assert(len);
    assert(grid);
    assert(grid->nv);
    assert(max_local_np > 0);
    assert(max_local_nm > 0);

    this->name = strdup(name);
    
    this->q = q;
    this->m = m;
    
    this->p = new PscParticle[max_local_np];
    this->max_np = max_local_np;
    
    this->pm = new PscParticleMover[max_local_nm];
    this->max_nm = max_local_nm;
    
    this->last_sorted       = INT64_MIN;
    this->sort_interval     = sort_interval;
    this->sort_out_of_place = sort_out_of_place;
    this->partition = new int[grid->nv + 1];
    
    this->grid_ = grid;
  }

  // ----------------------------------------------------------------------
  // PscSpecies::dtor
  
  ~PscSpecies()
  {
    delete[] partition;
    delete[] pm;
    delete[] p;
    free(name);
  }

  const Grid* grid() const { return grid_; }
  Grid* grid() { return grid_; }

  char * name;                        // Species name
  float q;                            // Species particle charge
  float m;                            // Species particle rest mass

  int np, max_np;                     // Number and max local particles
  PscParticle * ALIGNED(128) p;       // Array of particles for the species

  int nm, max_nm;                     // Number and max local movers in use
  PscParticleMover * ALIGNED(128) pm; // Particle movers

  int64_t last_sorted;                // Step when the particles were last
                                      // sorted.
  int sort_interval;                  // How often to sort the species
  int sort_out_of_place;              // Sort method
  int * ALIGNED(128) partition;       // Static array indexed 0:
  /**/                                // (nx+2)*(ny+2)*(nz+2).  Each value
  /**/                                // corresponds to the associated particle
  /**/                                // array index of the first particle in
  /**/                                // the cell.  Array is allocated and
  /**/                                // values computed in sort_p.  Purpose is
  /**/                                // for implementing collision models
  /**/                                // This is given in terms of the
  /**/                                // underlying's grids space filling
  /**/                                // curve indexing.  Thus, immediately
  /**/                                // after a sort:
  /**/                                //   sp->p[sp->partition[g->sfc[i]  ]:
  /**/                                //         sp->partition[g->sfc[i]+1]-1]
  /**/                                // are all the particles in voxel
  /**/                                // with local index i, while:
  /**/                                //   sp->p[ sp->partition[ j   ]:
  /**/                                //          sp->partition[ j+1 ] ]
  /**/                                // are all the particles in voxel
  /**/                                // with space filling curve index j.
  /**/                                // Note: SFC NOT IN USE RIGHT NOW THUS
  /**/                                // g->sfc[i]=i ABOVE.
private:
  Grid* grid_;                        // Underlying grid

public: // FIXME, shouldn't be public
  SpeciesId id;                       // Unique identifier for a species
  PscSpecies* next;                   // Next species in the list
};

// ======================================================================
// PscParticlesBase

template<class G, class BCL>
struct PscParticlesBase : public VpicListBase<PscSpecies<G>>
{
  typedef G Grid;
  typedef BCL ParticleBcList;
  typedef PscSpecies<Grid> Species;
  typedef VpicListBase<Species> Base;
  typedef PscParticle Particle;
  typedef PscParticleMover ParticleMover;

  using typename Base::iterator;
  using typename Base::const_iterator;
  using Base::size;
  using Base::begin;
  using Base::end;
  using Base::head_;
  using Base::push_front;

  static Species* create(const char * name, float q, float m,
			 int max_local_np, int max_local_nm,
			 int sort_interval, int sort_out_of_place,
			 Grid *grid)
  {
    return new Species(name, q, m, max_local_np, max_local_nm,
		       sort_interval, sort_out_of_place, grid);
  }

  size_t getNumSpecies() { return size(); }

  iterator find(int id)
  {
    return std::find_if(begin(), end(),
    			[&id](const Species &sp) { return sp.id == id; });
  }
  
  iterator find(const char *name)
  {
    return std::find_if(begin(), end(),
    			[&name](const Species &sp) { return strcmp(sp.name, name) == 0; });
  }
  
  Species *append(Species *sp)
  {
    assert(!sp->next);
    if (find(sp->name) != end()) {
      LOG_ERROR("There is already a species named \"%s\" in list", sp->name);
    }
    sp->id = size();
    push_front(*sp);
    return sp;
  }
  
  // ----------------------------------------------------------------------
  // inject_particle
  
  void inject_particle(PscParticlesBase& vmprts, const particle_inject& prt)
  {
    auto sp = vmprts.find(prt.kind);
    assert(sp != end());

    double x = prt.x[0], y = prt.x[1], z = prt.x[2];
    double ux = prt.u[0], uy = prt.u[1], uz = prt.u[2];
    double w = prt.w;
#if 0
    double age = 0.;
    int update_rhob = 0;
#endif

    int ix, iy, iz;

    const Grid* grid = sp->grid();
    const double x0 = (double)grid->x0, y0 = (double)grid->y0, z0 = (double)grid->z0;
    const double x1 = (double)grid->x1, y1 = (double)grid->y1, z1 = (double)grid->z1;
    const int    nx = grid->nx,         ny = grid->ny,         nz = grid->nz;
    
    // Do not inject if the particle is strictly outside the local domain
    // or if a far wall of local domain shared with a neighbor
    
    if ((x<x0) | (x>x1) | ( (x==x1) & (grid->bc[BOUNDARY(1,0,0)]>=0))) return;
    if ((y<y0) | (y>y1) | ( (y==y1) & (grid->bc[BOUNDARY(0,1,0)]>=0))) return;
    if ((z<z0) | (z>z1) | ( (z==z1) & (grid->bc[BOUNDARY(0,0,1)]>=0))) return;
    
    // This node should inject the particle
    
    if (sp->np>=sp->max_np) LOG_ERROR("No room to inject particle");
    
    // Compute the injection cell and coordinate in cell coordinate system
    // BJA:  Note the use of double precision here for accurate particle 
    //       placement on large meshes. 
 
    // The ifs allow for injection on the far walls of the local computational
    // domain when necessary
    
    x  = ((double)nx)*((x-x0)/(x1-x0)); // x is rigorously on [0,nx]
    ix = (int)x;                        // ix is rigorously on [0,nx]
    x -= (double)ix;                    // x is rigorously on [0,1)
    x  = (x+x)-1;                       // x is rigorously on [-1,1)
    if( ix==nx ) x = 1;                 // On far wall ... conditional move
    if( ix==nx ) ix = nx-1;             // On far wall ... conditional move
    ix++;                               // Adjust for mesh indexing

    y  = ((double)ny)*((y-y0)/(y1-y0)); // y is rigorously on [0,ny]
    iy = (int)y;                        // iy is rigorously on [0,ny]
    y -= (double)iy;                    // y is rigorously on [0,1)
    y  = (y+y)-1;                       // y is rigorously on [-1,1)
    if( iy==ny ) y = 1;                 // On far wall ... conditional move
    if( iy==ny ) iy = ny-1;             // On far wall ... conditional move
    iy++;                               // Adjust for mesh indexing

    z  = ((double)nz)*((z-z0)/(z1-z0)); // z is rigorously on [0,nz]
    iz = (int)z;                        // iz is rigorously on [0,nz]
    z -= (double)iz;                    // z is rigorously on [0,1)
    z  = (z+z)-1;                       // z is rigorously on [-1,1)
    if( iz==nz ) z = 1;                 // On far wall ... conditional move
    if( iz==nz ) iz = nz-1;             // On far wall ... conditional move
    iz++;                               // Adjust for mesh indexing

    Particle* p = sp->p + (sp->np++);
    p->dx = (float)x; // Note: Might be rounded to be on [-1,1]
    p->dy = (float)y; // Note: Might be rounded to be on [-1,1]
    p->dz = (float)z; // Note: Might be rounded to be on [-1,1]
    p->i  = VOXEL(ix,iy,iz, nx,ny,nz);
    p->ux = (float)ux;
    p->uy = (float)uy;
    p->uz = (float)uz;
    p->w  = w;

#if 0
    if (update_rhob) accumulate_rhob(fa, p, -sp->q);

    if (age!=0) {
      if( sp->nm >= sp->max_nm )
	LOG_WARN("No movers available to age injected  particle");
      ParticleMover * pm = sp->pm + sp->nm;
      age *= grid->cvac*grid->dt/sqrt( ux*ux + uy*uy + uz*uz + 1 );
      pm->dispx = ux*age*grid->rdx;
      pm->dispy = uy*age*grid->rdy;
      pm->dispz = uz*age*grid->rdz;
      pm->i     = sp->np-1;
      sp->nm += move_p( sp->p, pm, accumulator[0], grid, sp->q );
    }
#endif
  }
  
  Grid *grid()
  {
    assert(head_);
    return head_->grid();
  }
  
  Species* head()
  {
    return head_;
  }
};
 
#endif
