
#pragma once

template <typename Mparticles, typename MfieldsState,
          typename MfieldsInterpolator, typename MfieldsAccumulator,
          typename MfieldsHydro>
struct VpicParticlesOps
{
  using ParticleBcList = typename Mparticles::ParticleBcList;
  using FieldArray = typename MfieldsState::FieldArray;

  static void advance_p(Mparticles& mprts, MfieldsAccumulator& accumulator,
                        MfieldsInterpolator& interpolator)
  {
    auto prts = mprts[0];
    for (auto sp = prts.begin(); sp != prts.end(); ++sp) {
      TIC ::advance_p(&*sp, accumulator, interpolator.getPatch(0).ip());
      TOC(advance_p, 1);
    }
  }

  static void boundary_p(const ParticleBcList& pbc_list, Mparticles& mprts,
                         MfieldsState& mflds, MfieldsAccumulator& accumulator)
  {
    const particle_bc_t* pbc = pbc_list;
    ::boundary_p(const_cast<particle_bc_t*>(pbc), mprts.head(), mflds,
                 accumulator);
  }

  static void accumulate_rhob(MfieldsState& mflds, const particle_t* p,
                              float qsp)
  {
    FieldArray* fa = mflds;
    ::accumulate_rhob(fa->f, p, fa->g, qsp);
  }

  // ----------------------------------------------------------------------
  // drop_p

  static void drop_p(Mparticles& mprts, MfieldsState& mflds)
  {
    auto prts = mprts[0];
    for (auto sp = prts.begin(); sp != prts.end(); ++sp) {
      if (sp->nm) {
        LOG_WARN("Removing %i particles associated with unprocessed %s movers "
                 "(increase num_comm_round)",
                 sp->nm, sp->name);
      }
      // Drop the particles that have unprocessed movers due to a user defined
      // boundary condition. Particles of this type with unprocessed movers are
      // in the list of particles and move_p has set the voxel in the particle
      // to 8*voxel + face. This is an incorrect voxel index and in many cases
      // can in fact go out of bounds of the voxel indexing space. Removal is in
      // reverse order for back filling. Particle charge is accumulated to the
      // mesh before removing the particle.
      int nm = sp->nm;
      particle_mover_t* RESTRICT ALIGNED(16) pm = sp->pm + sp->nm - 1;
      particle_t* RESTRICT ALIGNED(128) p0 = sp->p;
      for (; nm; nm--, pm--) {
        int i = pm->i; // particle index we are removing
        p0[i].i >>= 3; // shift particle voxel down
        // accumulate the particle's charge to the mesh
        accumulate_rhob(mflds, p0 + i, sp->q);
        p0[i] = p0[sp->np - 1]; // put the last particle into position i
        sp->np--;               // decrement the number of particles
      }
      sp->nm = 0;
    }
  }

  static void uncenter_p(species_t* sp,
                         /*const*/ MfieldsInterpolator& interpolator)
  {
    ::uncenter_p(sp, interpolator.getPatch(0).ip());
  }

  static void sort_p(species_t* sp) { ::sort_p(sp); }

  static double energy_p(typename Mparticles::ConstSpeciesIterator sp,
                         const MfieldsInterpolator& interpolator)
  {
    return ::energy_p(&*sp, &interpolator);
  }
};
