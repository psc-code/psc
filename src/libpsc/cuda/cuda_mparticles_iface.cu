
#include "cuda_mparticles_iface.hxx"
#include "cuda_mparticles.cuh"
#include "bs.hxx"

template<typename BS>
cuda_mparticles<BS>* cuda_mparticles_iface<BS>::new_(const Grid_t& grid)
{
  return new cuda_mparticles(grid);
}

template<typename BS>
void cuda_mparticles_iface<BS>::delete_(cuda_mparticles* cmprts)
{
  delete cmprts;
}

template<typename BS>
int cuda_mparticles_iface<BS>::get_n_prts(cuda_mparticles* cmprts)
{
  return cmprts->get_n_prts();
}

template<typename BS>
std::vector<uint> cuda_mparticles_iface<BS>::get_size_all(const cuda_mparticles* cmprts)
{
  return cmprts->get_size_all();
}

template<typename BS>
void cuda_mparticles_iface<BS>::inject(cuda_mparticles* cmprts, const std::vector<cuda_mparticles_iface<BS>::Particle>& buf,
				       const std::vector<uint>& buf_n_by_patch)
{
  cmprts->inject(buf, buf_n_by_patch);
}

template<typename BS>
std::vector<uint> cuda_mparticles_iface<BS>::get_offsets(const cuda_mparticles* cmprts)
{
  return cmprts->get_offsets();
}

template<typename BS>
std::vector<typename cuda_mparticles_iface<BS>::Particle> cuda_mparticles_iface<BS>::get_particles(const cuda_mparticles* cmprts)
{
  return const_cast<cuda_mparticles*>(cmprts)->get_particles(); // FIXME cast
}

template<typename BS>
cuda_mparticles_iface<BS>::Particle cuda_mparticles_iface<BS>::get_particle(const cuda_mparticles* cmprts, int p, int n)
{
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


