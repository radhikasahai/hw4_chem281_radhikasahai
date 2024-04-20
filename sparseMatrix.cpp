#include <vector>
#include <iostream>
#include "vector.h"
#include "denseMatrix.h"

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
    SparseMatrix(int rows, int cols) : numRows(rows), numCols(cols) {
        data.resize(rows);
        diagonal.resize(rows, 0.0);
    }

    // // Function to insert non-zero element into the matrix
    // void insert(int i, int j, double value) {
    //     if (value != 0.0) {
    //         data[i].push_back({j, value}); // Insert non-zero element
    //         if (i == j) { //if same then diagonal 
    //             diagonal[i] = value; // Store diagonal elements separately
    //         }
    //     }
    // }
    void insert(int i, int j, double value) {
    if (value != 0.0) {
        RowElement element;
        element.j = j;
        element.d = value;
        data[i].push_back(element); // Insert non-zero element
        if (i == j) { //if same then diagonal 
            diagonal[i] = value; // Store diagonal elements separately
        }
    }
}

    // Function to get the number of rows
    int getNumRows() const {
        return numRows;
    }

    // Function to get the number of columns
    int getNumCols() const {
        return numCols;
    }

    // Function to access data in the matrix
    double operator()(int i, int j) const {
        if (i == j) {
            return diagonal[i]; // Return diagonal element
        } else {
            // Search for non-zero element in the row
            for (const RowElement& element : data[i]) {
                if (element.j == j) {
                    return element.d;
                }
            }
            return 0.0; // If no non-zero element found
        }
    }

    // Function to apply the sparse matrix to a vector
    std::vector<double> apply(const std::vector<double>& vector) const { // vector of dot product of input against each row 
        //we have to assume that "vector" is of length of num rows or cols

        std::vector<double> result(numRows);
        for (int i = 0; i < numRows; ++i) {
            double sum = 0.0;
            for (const RowElement& element : data[i]) {
                sum += element.d * vector[element.j];
            }
            sum += diagonal[i] * vector[i];
            result[i] = sum;
        }
        return result;
    }

    // Implicit conversion operator to dense matrix
    operator DenseMatrix() const {
        DenseMatrix dense(numRows, numCols);
        for (int i = 0; i < numRows; ++i) {
            for (const RowElement& element : data[i]) {
                dense(i, element.j) = element.d;
            }
            dense(i, i) = diagonal[i];
        }
        return dense;
    }
};


int main() {
    // Create a sparse matrix with 3 rows and 3 columns
    SparseMatrix sparse(3, 3);

    // Insert some non-zero elements
    sparse.insert(0, 0, 1.0); // Element at (0, 0)
    sparse.insert(1, 1, 2.0); // Element at (1, 1)
    sparse.insert(2, 2, 3.0); // Element at (2, 2)

    // Access and print the elements of the sparse matrix
    for (int i = 0; i < sparse.getNumRows(); ++i) {
        for (int j = 0; j < sparse.getNumCols(); ++j) {
            std::cout << sparse(i, j) << " ";
        }
        std::cout << std::endl;
    }

    // Convert the sparse matrix to a dense matrix
    DenseMatrix dense = sparse;

    // Print the dense matrix
    dense.print();

    return 0;
}
