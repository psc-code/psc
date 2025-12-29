#pragma once

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

/// @brief A parser that reads a dict-like map of parameter names to
/// parameter values from a single input file. The input file syntax is
/// extremely simple:
/// ```txt
/// param1 val1
/// param2 val2 thiswordisignored andsoisthis
/// ```
/// The first (space-separated) word in each line is interpreted as a parameter
/// name, and the second word is interpreted as a value. Words after the first
/// two are silently ignored, allowing for comments. If a parameter is
/// duplicated, all values but the last are silently ignored.
///
/// This class has no knowledge about what parameters should or shouldn't be
/// present, nor does it know what types its values should have. A value isn't
/// parsed to a specific type (e.g. `double`) until it is actually accessed as
/// that type by a user.
class ParsedParams
{
private:
  std::unordered_map<std::string, std::string> params;

public:
  ParsedParams(const std::string file_path)
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
      std::cout << "Failed to open params file: " << file_path << "\n";
      exit(EXIT_FAILURE);
    }
  }

  /// @brief Check if the parameter is present.
  /// @param paramName name of parameter
  /// @return whether or not the parameter is present
  bool has(const std::string paramName) { return params.count(paramName) == 1; }

  /// @brief Get a parameter, parsing it to the given type.
  /// @tparam T type of parameter
  /// @param paramName name of parameter
  /// @return the parameter
  template <typename T>
  T get(const std::string paramName);

  /// @brief Get a parameter if it's there, otherwise
  /// return the given default value.
  /// @tparam T type of parameter
  /// @param paramName name of parameter
  /// @param deflt default value of parameter
  /// @return the parameter
  template <typename T>
  T getOrDefault(const std::string paramName, T deflt)
  {
    if (has(paramName)) {
      return get<T>(paramName);
    }

    std::cout << "Warning: using default value for parameter '" << paramName
              << "': " << deflt << "\n";
    return deflt;
  }

  /// @brief Get a parameter and display a warning if it's there, otherwise
  /// return the given default value.
  /// @tparam T type of parameter
  /// @param paramName name of parameter
  /// @param deflt default value of parameter
  /// @return the parameter
  template <typename T>
  T getAndWarnOrDefault(const std::string paramName, T deflt)
  {
    if (!has(paramName)) {
      return deflt;
    }

    T val = get<T>(paramName);
    std::cout << "Warning: using non-default value for parameter '" << paramName
              << "': " << val << "\nDefault value of '" << deflt
              << "' is recommended.\n";
    return val;
  }

  /// @brief Display a warning if a parameter is present.
  /// @param paramName name of parameter
  /// @param advice user-friendly instructions on what to do instead
  /// @return whether or not the parameter was present
  bool warnIfPresent(const std::string paramName, const std::string advice)
  {
    if (!has(paramName)) {
      return false;
    }

    std::cout << "Warning: parameter " << paramName << " is deprecated.\n"
              << advice << "\n";
    return true;
  }

private:
  /// Retrieves an unparsed value, throwing a helpful error if the parameter is
  /// missing.
  std::string _getUnparsed(const std::string paramName)
  {
    if (has(paramName)) {
      return params.at(paramName);
    }

    throw std::out_of_range("missing required input parameter: " + paramName);
  }
};

// get implementations

template <>
bool ParsedParams::get<bool>(const std::string paramName)
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
double ParsedParams::get<double>(const std::string paramName)
{
  return std::stod(_getUnparsed(paramName));
}

template <>
int ParsedParams::get<int>(const std::string paramName)
{
  return std::stoi(_getUnparsed(paramName));
}

template <>
float ParsedParams::get<float>(const std::string paramName)
{
  return std::stof(_getUnparsed(paramName));
}

template <>
std::string ParsedParams::get<std::string>(const std::string paramName)
{
  return _getUnparsed(paramName);
}
