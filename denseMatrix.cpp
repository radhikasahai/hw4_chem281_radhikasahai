// #include <vector>
// #include <iostream>
// #include "denseMatrix.h"

// class DenseMatrix {
// private:
//     std::vector<std::vector<double> > data;
//     int numRows;
//     int numCols;

// public:
//     // Constructors
//     DenseMatrix() : numRows(0), numCols(0) {}

//     DenseMatrix(int rows, int cols, double value = 0.0) : numRows(rows), numCols(cols) {
//         data.resize(rows, std::vector<double>(cols, value));
//     }

//     // Destructor
//     ~DenseMatrix() {}

//     // Function to get the number of rows
//     int getNumRows() const {
//         return numRows;
//     }

//     // Function to get the number of columns
//     int getNumCols() const {
//         return numCols;
//     }

//     // Function to access elements by indices
//     double& operator()(int row, int col) {
//         return data[row][col];
//     }

//     // Function to access elements by indices (const version)
//     const double& operator()(int row, int col) const {
//         return data[row][col];
//     }

//     // Function to resize the matrix
//     void resize(int rows, int cols, double value = 0.0) {
//         data.resize(rows, std::vector<double>(cols, value));
//         numRows = rows;
//         numCols = cols;
//     }

//     // Function to clear the matrix
//     void clear() {
//         data.clear();
//         numRows = 0;
//         numCols = 0;
//     }

//     // Function to print the matrix
//     void print() const {
//         for (int i = 0; i < numRows; ++i) {
//             for (int j = 0; j < numCols; ++j) {
//                 std::cout << data[i][j] << " ";
//             }
//             std::cout << std::endl;
//         }
//     }
// };

#include "denseMatrix.h"

// Constructors
DenseMatrix::DenseMatrix() : numRows(0), numCols(0) {}

DenseMatrix::DenseMatrix(int rows, int cols, double value) : numRows(rows), numCols(cols) {
    data.resize(rows, std::vector<double>(cols, value));
}

// Destructor
DenseMatrix::~DenseMatrix() {}

// Function to get the number of rows
int DenseMatrix::getNumRows() const {
    return numRows;
}

// Function to get the number of columns
int DenseMatrix::getNumCols() const {
    return numCols;
}

// Function to access elements by indices
double& DenseMatrix::operator()(int row, int col) {
    return data[row][col];
}

const double& DenseMatrix::operator()(int row, int col) const {
    return data[row][col];
}

// Function to resize the matrix
void DenseMatrix::resize(int rows, int cols, double value) {
    data.resize(rows, std::vector<double>(cols, value));
    numRows = rows;
    numCols = cols;
}

// Function to clear the matrix
void DenseMatrix::clear() {
    data.clear();
    numRows = 0;
    numCols = 0;
}

// Function to print the matrix
void DenseMatrix::print() const {
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


// int main() {
//     // Create a 5x5 matrix of 0s
//     DenseMatrix matrix(5, 5);

//     // Print the matrix
//     std::cout << "Matrix of zeros:" << std::endl;
//     matrix.print();

//     return 0;
// }