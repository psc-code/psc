
#include "cuda_mparticles_iface.hxx"
#include "cuda_mparticles.cuh"
#include "bs.hxx"

#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

// ======================================================================
// cuda_mparticles_iface implementation

template<typename BS>
cuda_mparticles<BS>* cuda_mparticles_iface<BS>::new_(const Grid_t& grid)
{
  dprintf("CMPRTS: new\n");
  return new cuda_mparticles(grid);
}

template<typename BS>
void cuda_mparticles_iface<BS>::delete_(cuda_mparticles* cmprts)
{
  dprintf("CMPRTS: delete\n");
  delete cmprts;
}

template<typename BS>
int cuda_mparticles_iface<BS>::get_n_prts(cuda_mparticles* cmprts)
{
  dprintf("CMPRTS: get_n_prts\n");
  return cmprts->get_n_prts();
}

template<typename BS>
std::vector<uint> cuda_mparticles_iface<BS>::get_size_all(const cuda_mparticles* cmprts)
{
  dprintf("CMPRTS: get_size_all\n");
  auto n_prts_by_patch = cmprts->get_size_all();
  for (int p = 0; p < n_prts_by_patch.size(); p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
  return n_prts_by_patch;
}

template<typename BS>
void cuda_mparticles_iface<BS>::inject(cuda_mparticles* cmprts, const std::vector<cuda_mparticles_iface<BS>::Particle>& buf,
				       const std::vector<uint>& buf_n_by_patch)
{
  dprintf("CMPRTS: inject\n");
  cmprts->inject(buf, buf_n_by_patch);
}

template<typename BS>
std::vector<uint> cuda_mparticles_iface<BS>::get_offsets(const cuda_mparticles* cmprts)
{
  dprintf("CMPRTS: get_offsets\n");
  return cmprts->get_offsets();
}

template<typename BS>
std::vector<typename cuda_mparticles_iface<BS>::Particle> cuda_mparticles_iface<BS>::get_particles(const cuda_mparticles* cmprts)
{
  dprintf("CMPRTS: get_particles\n");
  return const_cast<cuda_mparticles*>(cmprts)->get_particles(); // FIXME cast
}

template<typename BS>
cuda_mparticles_iface<BS>::Particle cuda_mparticles_iface<BS>::get_particle(const cuda_mparticles* cmprts, int p, int n)
{
  dprintf("CMPRTS: get_particle\n");
  return const_cast<cuda_mparticles*>(cmprts)->get_particle(p, n); // FIXME cast
}

template<typename BS>
bool cuda_mparticles_iface<BS>::check_after_push(cuda_mparticles* cmprts)
{
  return cmprts->check_bidx_after_push(); // FIXME, naming inconsistence
}

template<typename BS>
void cuda_mparticles_iface<BS>::dump(const cuda_mparticles* cmprts, const std::string& filename)
{
  cmprts->dump(filename);
}

template struct cuda_mparticles_iface<BS144>;
template struct cuda_mparticles_iface<BS444>;


