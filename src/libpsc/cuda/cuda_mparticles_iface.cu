
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
void cuda_mparticles_iface<BS>::inject(cuda_mparticles* cmprts, const std::vector<ParticleSimple<real_t>>& buf,
				       const std::vector<uint>& buf_n_by_patch)
{
  cmprts->inject(buf, buf_n_by_patch);
}

template struct cuda_mparticles_iface<BS144>;
template struct cuda_mparticles_iface<BS444>;


