#pragma once

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <mpi.h>

#include "../libpsc/vpic/psc_vpic_bits.h"

namespace
{
template <typename T>
std::string to_str(const T& val)
{
  std::stringstream ss;
  ss << val;
  return ss.str();
}

int rank()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}
} // namespace

/**
 * @brief A parser that reads a dict-like map of parameter names to
 * parameter values from a single input file. The input file syntax is
 * extremely simple:
 * ```txt
 * param1 val1
 * param2 val2 thiswordisignored andsoisthis
 * ```
 * The first (space-separated) word in each line is interpreted as a parameter
 * name, and the second word is interpreted as a value. Words after the first
 * two are silently ignored, allowing for comments. If a parameter is
 * duplicated, all values but the last are silently ignored.
 *
 * This class has no knowledge about what parameters should or shouldn't be
 * present, nor does it know what types its values should have. A value isn't
 * parsed to a specific type (e.g. `double`) until it is actually accessed as
 * that type by a user.
 */
class InputParams
{
private:
  std::unordered_map<std::string, std::string> params;

public:
  InputParams(const std::string file_path)
  {
    // iterate over each line
    std::ifstream ifs(file_path);

    if (ifs.is_open()) {
      for (std::string line; std::getline(ifs, line);) {

        // parse first two words within line
        std::istringstream iss(line);
        std::string paramName, paramVal;
        if (iss >> paramName >> paramVal)
          params[paramName] = paramVal;
      }

      ifs.close();
    } else {
      LOG_ERROR("Failed to open parameter file: %s\n", file_path.c_str());
    }
  }

  /**
   * @brief Check if the parameter is present.
   * @param paramName paramName name of parameter
   * @return whether or not the parameter is present
   */
  bool has(const std::string paramName) { return params.count(paramName) == 1; }

  /**
   * @brief Get a parameter, parsing it to the given type.
   * @tparam T type of parameter
   * @param paramName name of parameter
   * @return the parameter
   */
  template <typename T>
  T get(const std::string paramName)
  {
    try {
      return _getParsed<T>(paramName);
    } catch (const std::invalid_argument& e) {
      std::string unparsed = _getUnparsed(paramName);
      // TODO ensure human-readable type name
      LOG_ERROR(
        "Unable to parse parameter '%s', which has value '%s', to type %s\n",
        paramName.c_str(), unparsed.c_str(), typeid(T).name());
    }
  }

  /**
   * @brief Get a parameter if it's there, otherwise return the given default
   * value.
   * @tparam T type of parameter
   * @param paramName name of parameter
   * @param deflt default value of parameter
   * @return the parameter
   */
  template <typename T>
  T getOrDefault(const std::string paramName, T deflt)
  {
    if (has(paramName)) {
      return get<T>(paramName);
    }

    if (rank() == 0) {
      LOG_WARN("Using default value for parameter '%s': %s\n",
               paramName.c_str(), to_str(deflt).c_str());
    }

    return deflt;
  }

  /**
   * @brief Get a parameter and display a warning if it's there, otherwise
   * return the given default value.
   * @tparam T type of parameter
   * @param paramName name of parameter
   * @param deflt default value of parameter
   * @return the parameter
   */
  template <typename T>
  T getAndWarnOrDefault(const std::string paramName, T deflt)
  {
    if (!has(paramName)) {
      return deflt;
    }

    T val = get<T>(paramName);

    if (rank() == 0) {
      LOG_WARN(
        "Using non-default value for parameter '%s': %s (default value of "
        "%s is recommended)\n",
        paramName.c_str(), t_c_str(val), to_str(deflt).c_str());
    }

    return val;
  }

  /**
   * @brief Display a warning if a parameter is present.
   * @param paramName name of parameter
   * @param advice user-friendly instructions on what to do instead
   * @return whether or not the parameter was present
   */
  bool warnIfPresent(const std::string paramName, const std::string advice)
  {
    if (!has(paramName)) {
      return false;
    }

    if (rank() == 0) {
      LOG_WARN("Parameter '%s' is deprecated; %s\n", paramName.c_str(),
               advice.c_str());
    }

    return true;
  }

private:
  /**
   * @brief Retrieves an unparsed value, throwing a helpful error if the
   * parameter is not present.
   * @param paramName name of parameter
   * @return the unparsed value (as a string)
   */
  std::string _getUnparsed(const std::string paramName)
  {
    if (has(paramName)) {
      return params.at(paramName);
    }

    LOG_ERROR("Required parameter '%s' is absent\n", paramName.c_str());
  }

  /**
   * @brief Retrieve and parse a value, possibly throwing
   * `std::invalid_argument`.
   * @tparam T the type of the parameter value
   * @param paramName name of parameter
   * @return the parsed value
   */
  template <typename T>
  T _getParsed(const std::string paramName);
};

// get implementations

template <>
bool InputParams::_getParsed<bool>(const std::string paramName)
{
  auto lowercase = _getUnparsed(paramName);
  std::transform(lowercase.begin(), lowercase.end(), lowercase.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (lowercase == "true") {
    return true;
  } else if (lowercase == "false") {
    return false;
  } else {
    throw std::invalid_argument(_getUnparsed(paramName));
  }
}

template <>
double InputParams::_getParsed<double>(const std::string paramName)
{
  return std::stod(_getUnparsed(paramName));
}

template <>
int InputParams::_getParsed<int>(const std::string paramName)
{
  return std::stoi(_getUnparsed(paramName));
}

template <>
float InputParams::_getParsed<float>(const std::string paramName)
{
  return std::stof(_getUnparsed(paramName));
}

template <>
std::string InputParams::_getParsed<std::string>(const std::string paramName)
{
  return _getUnparsed(paramName);
}
