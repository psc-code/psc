
#pragma once

#include <mrc_common.h>

#include <deque>
#include <iostream>

#include "File.h"

namespace kg
{
namespace io
{

template <typename T>
class Variable;

// ======================================================================
// Engine

class Engine
{
public:
  Engine(File&& file, MPI_Comm comm);

  // ----------------------------------------------------------------------
  // put

  template <class T, class... Args>
  void put(const std::string& pfx, const T& datum, Args&&... args);

  template <class T>
  void putLocal(const std::string& pfx, const T& datum,
                Mode launch = Mode::NonBlocking);

  template <template <typename...> class Var, class T, class... Args>
  void put(const std::string& pfx, const T& datum, Args&&... args);

  // ----------------------------------------------------------------------
  // get

  template <class T, class... Args>
  void get(const std::string& pfx, T& datum, Args&&... args);

  template <class T>
  void getLocal(const std::string& pfx, T& datum, Mode launch = Mode::NonBlocking);

  template <template <typename...> class Var, class T, class... Args>
  void get(const std::string& pfx, T& datum, Args&&... args);

  // ----------------------------------------------------------------------
  // performPuts

  void performPuts();

  // ----------------------------------------------------------------------
  // performGets

  void performGets();

  // ----------------------------------------------------------------------
  // variableShape

  template <typename T>
  Dims variableShape();

  // ----------------------------------------------------------------------
  // internal

  template <typename T>
  void putVariable(const T* data, const Mode launch, const Dims& shape,
                   const Extents& selection = {},
                   const Extents& memory_selection = {});

  template <typename T>
  void getVariable(T* data, const Mode launch, const Extents& selection = {},
                   const Extents& memory_selection = {});

  template <typename T>
  void putAttribute(const T& datum);

  template <typename T>
  void putAttribute(const T* data, size_t size);

  template <typename T>
  void getAttribute(T& datum);

  template <typename T>
  void getAttribute(std::vector<T>& data);

  // ----------------------------------------------------------------------
  // close

  void close();

  int mpiRank() const;
  int mpiSize() const;

  std::string prefix() const
  {
    std::string s;
    bool first = true;
    for (auto& pfx : prefixes_) {
      if (!first) {
        s += "::";
      }
      s += pfx;
      first = false;
    }
    return s;
  }

private:
  File file_;
  std::deque<std::string> prefixes_;
  int mpi_rank_;
  int mpi_size_;
};

} // namespace io
} // namespace kg

#include "Engine.inl"
