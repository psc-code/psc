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

template <int n_rows, int n_cols>
class Parsed
{
private:
  double data[n_rows][n_cols];
  const int indep_col;
  double indep_var_step;

  // ----------------------------------------------------------------------
  // get_row
  // Gets index of row containing greatest lower bound of given indep_var_val

  int get_row(double indep_var_val)
  {
    // initial guess; should be precise assuming indep_var is linearly spaced
    int row = std::max(int(indep_var_val / indep_var_step), n_rows - 1);
    while (indep_var_val < data[row][indep_col])
      row--;
    return row;
  }

public:
  // ----------------------------------------------------------------------
  // ctor
  // Parses all the data

  Parsed(int indep_col) : indep_col(indep_col) {}

  void loadData(const std::string file_path, int lines_to_skip)
  {
    std::ifstream file(file_path);

    for (int i = 0; i < lines_to_skip; i++)
      file.ignore(512, '\n');

    // iterate over each line
    int row = 0;

    for (std::string line; std::getline(file, line);) {
      assert(row < n_rows);

      // iterate over each entry within a line
      std::istringstream iss(line);
      int col = 0;

      for (std::string result; iss >> result;) {
        assert(col < n_cols);

        // write entry to data
        data[row][col] = std::stod(result);
        col++;
      }
      assert(col == n_cols);
      row++;
    }
    assert(row == n_rows);
    file.close();

    indep_var_step = data[1][indep_col] - data[0][indep_col];
  }

  // ----------------------------------------------------------------------
  // get_interpolated
  // Calculates and returns a linearly interpolated value (specified by col) at
  // given value of independent variable

  double get_interpolated(int col, double indep_var_val)
  {
    assert(indep_var_val >= data[0][indep_col]);
    assert(indep_var_val <= data[n_rows - 1][indep_col]);

    int row = get_row(indep_var_val);

    if (data[row][indep_col] == indep_var_val)
      return data[row][col];

    // weights for linear interpolation
    double w1 = data[row + 1][indep_col] - indep_var_val;
    double w2 = indep_var_val - data[row][indep_col];

    return (w1 * data[row][col] + w2 * data[row + 1][col]) / (w1 + w2);
  }
};