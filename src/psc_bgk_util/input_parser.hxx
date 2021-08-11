#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

// ======================================================================
// Parsed
// Parses a space-separated list of values, such as a tsv file.
// Assuming there is a single independent variable, linearly interpolates other
// values.

class Parsed
{
private:
  const int nrows, ncols;
  std::vector<double> data;
  const int indep_col;
  double indep_var_step;

  // ----------------------------------------------------------------------
  // get_row
  // Gets index of row containing greatest lower bound of given indep_var_val

  int get_row(double indep_var_val)
  {
    // initial guess; should be precise assuming indep_var is linearly spaced
    int row = std::max((int)(indep_var_val / indep_var_step), nrows - 1);
    while (indep_var_val < (*this)[row][indep_col])
      row--;
    return row;
  }

public:
  // ----------------------------------------------------------------------
  // ctor

  Parsed(int nrows, int ncols, int indep_col)
    : nrows(nrows), ncols(ncols), data(nrows * ncols), indep_col(indep_col)
  {}

  // ----------------------------------------------------------------------
  // loadData
  // Parses all the data

  void loadData(const std::string file_path, int lines_to_skip)
  {
    std::ifstream file(file_path);

    if (!file.is_open()) {
      std::cout << "Failed to open input file: " << file_path << std::endl;
      exit(EXIT_FAILURE);
    }

    for (int i = 0; i < lines_to_skip; i++)
      file.ignore(512, '\n');

    // iterate over each line
    int row = 0;

    for (std::string line; std::getline(file, line);) {
      assert(row < nrows);

      // iterate over each entry within a line
      std::istringstream iss(line);
      int col = 0;

      for (std::string result; iss >> result;) {
        assert(col < ncols);

        // write entry to data
        (*this)[row][col] = std::stod(result);
        col++;
      }
      assert(col == ncols);
      row++;
    }
    assert(row == nrows);
    file.close();

    indep_var_step = (*this)[1][indep_col] - (*this)[0][indep_col];
  }

  // ----------------------------------------------------------------------
  // get_interpolated
  // Calculates and returns a linearly interpolated value (specified by col) at
  // given value of independent variable

  double get_interpolated(int col, double indep_var_val)
  {
    assert(indep_var_val >= (*this)[0][indep_col]);
    assert(indep_var_val <= (*this)[nrows - 1][indep_col]);

    int row = get_row(indep_var_val);

    if ((*this)[row][indep_col] == indep_var_val)
      return (*this)[row][col];

    // weights for linear interpolation
    double w1 = (*this)[row + 1][indep_col] - indep_var_val;
    double w2 = indep_var_val - (*this)[row][indep_col];

    return (w1 * (*this)[row][col] + w2 * (*this)[row + 1][col]) / (w1 + w2);
  }

  double* operator[](const int row) { return data.data() + row * ncols; }
};