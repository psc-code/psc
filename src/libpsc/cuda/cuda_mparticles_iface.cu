
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
cuda_mparticles_iface<BS>::CudaMparticles* cuda_mparticles_iface<BS>::new_(const Grid_t& grid)
{
  dprintf("CMPRTS: new\n");
  return new CudaMparticles(grid);
}

template<typename BS>
void cuda_mparticles_iface<BS>::delete_(CudaMparticles* cmprts)
{
  dprintf("CMPRTS: delete\n");
  delete cmprts;
}

template<typename BS>
int cuda_mparticles_iface<BS>::size(CudaMparticles* cmprts)
{
  dprintf("CMPRTS: size\n");
  return cmprts->size();
}

template<typename BS>
std::vector<uint> cuda_mparticles_iface<BS>::sizeByPatch(const CudaMparticles* cmprts)
{
  dprintf("CMPRTS: sizeByPatch\n");
  auto n_prts_by_patch = cmprts->sizeByPatch();
  for (int p = 0; p < n_prts_by_patch.size(); p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
  return n_prts_by_patch;
}

template<typename BS>
void cuda_mparticles_iface<BS>::inject(CudaMparticles* cmprts, const std::vector<cuda_mparticles_iface<BS>::Particle>& buf,
				       const std::vector<uint>& buf_n_by_patch)
{
  dprintf("CMPRTS: inject\n");
  cmprts->inject(buf, buf_n_by_patch);
}

template<typename BS>
std::vector<uint> cuda_mparticles_iface<BS>::get_offsets(const CudaMparticles* cmprts)
{
  dprintf("CMPRTS: get_offsets\n");
  return cmprts->get_offsets();
}

template<typename BS>
std::vector<typename cuda_mparticles_iface<BS>::Particle> cuda_mparticles_iface<BS>::get_particles(const CudaMparticles* cmprts)
{
  dprintf("CMPRTS: get_particles\n");
  return const_cast<CudaMparticles*>(cmprts)->get_particles(); // FIXME cast
}

template<typename BS>
cuda_mparticles_iface<BS>::Particle cuda_mparticles_iface<BS>::get_particle(const CudaMparticles* cmprts, int p, int n)
{
  dprintf("CMPRTS: get_particle\n");
  return const_cast<CudaMparticles*>(cmprts)->get_particle(p, n); // FIXME cast
}

template<typename BS>
bool cuda_mparticles_iface<BS>::check_after_push(CudaMparticles* cmprts)
{
  return cmprts->check_bidx_after_push(); // FIXME, naming inconsistence
}
template<typename BS>
bool cuda_mparticles_iface<BS>::need_reorder(CudaMparticles* cmprts)
{
  return cmprts->need_reorder;
}
template<typename BS>
void cuda_mparticles_iface<BS>::dump(const CudaMparticles* cmprts, const std::string& filename)
{
  cmprts->dump(filename);
}

template struct cuda_mparticles_iface<BS144>;
template struct cuda_mparticles_iface<BS444>;


