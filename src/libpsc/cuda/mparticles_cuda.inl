
#include <kg/io.h>

template <typename BS>
class kg::io::Descr<MparticlesCuda<BS>>
{
public:
  using Mparticles = MparticlesCuda<BS>;

  void put(kg::io::Engine& writer, const Mparticles& mprts_cuda_,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    auto& mprts_cuda = const_cast<Mparticles&>(mprts_cuda_); // FIXME
    auto& mprts = mprts_cuda.template get_as<MparticlesSingle>();

    writer.put("mprts", mprts, launch);
    writer.performPuts();

    mprts_cuda.put_as(mprts, MP_DONT_COPY);
  }

  void get(kg::io::Engine& reader, Mparticles& mprts_cuda,
           const kg::io::Mode launch = kg::io::Mode::NonBlocking)
  {
    auto& mprts = mprts_cuda.template get_as<MparticlesSingle>();

    reader.get("mprts", mprts, launch);
    reader.performGets();

    mprts_cuda.put_as(mprts);
  }
};
