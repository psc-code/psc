#pragma once

#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

class ParsedParams
{
private:
  std::unordered_map<std::string, std::string> params;

public:
  ParsedParams(const std::string file_path)
  {
    // iterate over each line
    std::ifstream file(file_path);
    for (std::string line; std::getline(file, line);) {

      // parse first two words within line
      std::istringstream iss(line);
      std::string paramName, paramVal;
      std::getline(iss, paramName, ' ');
      std::getline(iss, paramVal, ' ');

      params[paramName] = paramVal;
    }

    file.close();
  }

  template <typename T>
  T get(const std::string paramName);
};

// implementations

template <>
bool ParsedParams::get<bool>(const std::string paramName)
{
  bool b;
  std::istringstream(params[paramName]) >> std::boolalpha >> b;
  return b;
}

template <>
double ParsedParams::get<double>(const std::string paramName)
{
  return std::stod(params[paramName]);
}

template <>
int ParsedParams::get<int>(const std::string paramName)
{
  return std::stoi(params[paramName]);
}

template <>
float ParsedParams::get<float>(const std::string paramName)
{
  return std::stof(params[paramName]);
}

template <>
std::string ParsedParams::get<std::string>(const std::string paramName)
{
  return params[paramName];
}