#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include "denseMatrix.h" // Include the header file for DenseMatrix

// Struct to represent non-zero elements in a row
struct RowElement {
    int j; // Column index
    double d; // Value
};

class SparseMatrix {
private:
    std::vector<std::vector<RowElement> > data; // Underlying data structure
    std::vector<double> diagonal; // Dense vector for diagonal elements
    int numRows;
    int numCols;

public:
    // Constructors
    SparseMatrix(int rows, int cols);

    // Function to insert non-zero element into the matrix
    void insert(int i, int j, double value);

    // Function to get the number of rows
    int getNumRows() const;

    // Function to get the number of columns
    int getNumCols() const;

    // Function to access data in the matrix
    double operator()(int i, int j) const;

    // Function to apply the sparse matrix to a vector
    std::vector<double> apply(const std::vector<double>& vector) const;

    // Implicit conversion operator to dense matrix
    operator DenseMatrix() const;
};

#endif // SPARSE_MATRIX_H

