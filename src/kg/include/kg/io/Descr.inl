
#include "Engine.h"

namespace kg
{
namespace io
{

// ======================================================================
// Descr<T>

template <class T, class Enable>
inline void Descr<T, Enable>::put(Engine& writer, const T& value, Mode launch)
{
  writer.putAttribute(value);
}

template <class T, class Enable>
inline void Descr<T, Enable>::get(Engine& reader, T& value, Mode launch)
{
  reader.getAttribute(value);
}

// ======================================================================
// Descr<T[N]>

template <class T, size_t N>
inline void Descr<T[N]>::put(Engine& writer, const Descr<T[N]>::value_type& arr,
                             Mode launch)
{
  writer.putAttribute(arr, N);
}

template <class T, size_t N>
inline void Descr<T[N]>::get(Engine& reader, Descr<T[N]>::value_type& arr,
                             Mode launch)
{
  std::vector<T> vals;
  reader.getAttribute(vals);
  assert(vals.size() == N);
  std::copy(vals.begin(), vals.end(), arr);
}

// ======================================================================
// Descr<std::vector<T>>

template <class T>
inline void Descr<std::vector<T>>::put(Engine& writer,
                                       const std::vector<T>& vec, Mode launch)
{
  writer.putAttribute(vec.data(), vec.size());
}

template <class T>
inline void Descr<std::vector<T>>::get(Engine& reader, std::vector<T>& vec,
                                       Mode launch)
{
  reader.getAttribute(vec);
}

// ======================================================================
// Descr<Vec3<T>>

template <class T>
inline void Descr<Vec3<T>>::put(Engine& writer, const Vec3<T>& vec, Mode launch)
{
  writer.putAttribute(vec.data(), 3);
}

template <class T>
inline void Descr<Vec3<T>>::get(Engine& reader, Vec3<T>& vec, Mode launch)
{
  std::vector<T> vals;
  reader.getAttribute(vals);
  assert(vals.size() == 3);
  vec = {vals[0], vals[1], vals[2]};
}

// ======================================================================
// Local<T>

template <class T, class Enable>
inline void Local<T, Enable>::put(Engine& writer, const T& value, Mode launch)
{
  writer.putVariable(&value, launch, {LocalValueDim});
}

template <class T, class Enable>
inline void Local<T, Enable>::get(Engine& reader, T& value, Mode launch)
{
  auto shape = reader.variableShape<T>();
  assert(shape == Dims{static_cast<size_t>(reader.mpiSize())});

  // FIXME, setSelection doesn't work, so read the whole thing
  std::vector<T> vals(shape[0]);
  reader.getVariable(vals.data(), launch);
  value = vals[reader.mpiRank()];
}

} // namespace io
} // namespace kg
