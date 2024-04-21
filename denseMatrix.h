#pragma once

#include <vector>
#include <iostream>

class DenseMatrix {
private:
    std::vector<std::vector<double> > data;
    int numRows;
    int numCols;

public:
    // Constructors
    DenseMatrix();
    DenseMatrix(int rows, int cols, double value = 0.0);

    // Destructor
    ~DenseMatrix();

    // Function to get the number of rows
    int getNumRows() const;

    // Function to get the number of columns
    int getNumCols() const;

    // Function to access elements by indices
    double& operator()(int row, int col);
    const double& operator()(int row, int col) const;

    // Function to resize the matrix
    void resize(int rows, int cols, double value = 0.0);

    // Function to clear the matrix
    void clear();

    // Function to print the matrix
    void print() const;
};

