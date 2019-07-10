
#include "Descr.h"

namespace kg
{
namespace io
{

// ======================================================================
// Engine

inline Engine::Engine(File&& file, MPI_Comm comm)
  : file_{std::move(file)}
{
  MPI_Comm_rank(comm, &mpi_rank_);
  MPI_Comm_size(comm, &mpi_size_);
}

// ----------------------------------------------------------------------
// put

template <class T, class... Args>
inline void Engine::put(const std::string& pfx, const T& datum, Args&&... args)
{
  put<Descr>(pfx, datum, std::forward<Args>(args)...);
}

template <class T>
inline void Engine::putLocal(const std::string& pfx, const T& datum,
                             Mode launch)
{
  put<Local>(pfx, datum, launch);
}

template <template <typename...> class Var, class T, class... Args>
inline void Engine::put(const std::string& pfx, const T& datum, Args&&... args)
{
  prefixes_.push_back(pfx);
  Var<T> var;
  var.put(*this, datum, std::forward<Args>(args)...);
  prefixes_.pop_back();
}

// ----------------------------------------------------------------------
// get

template <class T>
inline void Engine::getLocal(const std::string& pfx, T& datum, Mode launch)
{
  get<Local>(pfx, datum, launch);
}

template <class T, class... Args>
inline void Engine::get(const std::string& pfx, T& datum, Args&&... args)
{
  get<Descr>(pfx, datum, std::forward<Args>(args)...);
}

template <template <typename...> class Var, class T, class... Args>
inline void Engine::get(const std::string& pfx, T& datum, Args&&... args)
{
  prefixes_.push_back(pfx);
  // mprintf("get<Var> pfx %s -- %s\n", pfx.c_str(), prefix().c_str());
  Var<T> var;
  var.get(*this, datum, std::forward<Args>(args)...);
  prefixes_.pop_back();
}

// ----------------------------------------------------------------------
// performPuts

inline void Engine::performPuts()
{
  file_.performPuts();
}

// ----------------------------------------------------------------------
// performGets

inline void Engine::performGets()
{
  file_.performGets();
}

// ----------------------------------------------------------------------
// internal

template <typename T>
inline void Engine::putVariable(const T* data, const Mode launch,
                                const Dims& shape, const Extents& selection,
                                const Extents& memory_selection)
{
  file_.putVariable(prefix(), data, launch, shape, selection, memory_selection);
}

template <typename T>
inline void Engine::putAttribute(const T& datum)
{
  file_.putAttribute(prefix(), &datum, 1);
}

template <typename T>
inline void Engine::putAttribute(const T* data, size_t size)
{
  file_.putAttribute(prefix(), data, size);
}

template <typename T>
inline void Engine::getVariable(T* data, const Mode launch,
                                const Extents& selection,
                                const Extents& memory_selection)
{
  file_.getVariable(prefix(), data, launch, selection, memory_selection);
}

template <typename T>
inline void Engine::getAttribute(T& datum)
{
  auto size = file_.sizeAttribute(prefix());
  assert(size == 1);
  file_.getAttribute(prefix(), &datum);
}

template <typename T>
inline void Engine::getAttribute(std::vector<T>& vec)
{
  auto size = file_.sizeAttribute(prefix());
  vec.resize(size);
  file_.getAttribute(prefix(), vec.data());
}

// ----------------------------------------------------------------------
// variableShape

template <typename T>
inline Dims Engine::variableShape()
{
  return file_.shapeVariable(prefix());
}

// ----------------------------------------------------------------------
// close

inline void Engine::close()
{
  file_.close();
}

inline int Engine::mpiRank() const
{
  return mpi_rank_;
}

inline int Engine::mpiSize() const
{
  return mpi_size_;
}

} // namespace io
} // namespace kg
