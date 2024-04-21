#include <iostream>
#include <vector>
#include <tuple>
#include <lapacke.h>  
// g++ -o lapack lapack.cpp -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas


class Vector { //Vector object w/ pointer access for lapack
        std::vector<double> data;

    public:
        Vector(int n) : data(n) {}

        double* access_raw_pointer() { return data.data(); }
        int size() const { return data.size(); }
};

class DenseMatrix {
private:
    std::vector<std::vector<double> > data;
    int numRows;
    int numCols;

public:
    // Constructors
    DenseMatrix() : numRows(0), numCols(0) {}

    DenseMatrix(int rows, int cols, double value = 0.0) : numRows(rows), numCols(cols) {
        data.resize(rows, std::vector<double>(cols, value));
    }

    // Destructor
    ~DenseMatrix() {}

    // Function to get the number of rows
    int getNumRows() const {
        return numRows;
    }

    // Function to get the number of columns
    int getNumCols() const {
        return numCols;
    }

    // double* access_raw_pointer() { return data.data(); }

    // Function to get the underlying storage of a vector of vectors
    double* access_raw_pointer() {
        // Check if the outer vector is not empty
        if (data.empty()) {
            return nullptr; // Return nullptr if the outer vector is empty
        }

        // Access the first inner vector
        std::vector<double>& firstInnerVec = data[0];

        // Check if the inner vector is not empty
        if (firstInnerVec.empty()) {
            return nullptr; // Return nullptr if the inner vector is empty
        }

        // Return a pointer to the beginning of the underlying storage of the first inner vector
        return &firstInnerVec[0];
    }

    // Function to access elements by indices
    double& operator()(int row, int col) {
        return data[row][col];
    }

    // Function to access elements by indices (const version)
    const double& operator()(int row, int col) const {
        return data[row][col];
    }

    // Function to resize the matrix
    void resize(int rows, int cols, double value = 0.0) {
        data.resize(rows, std::vector<double>(cols, value));
        numRows = rows;
        numCols = cols;
    }

    // Function to clear the matrix
    void clear() {
        data.clear();
        numRows = 0;
        numCols = 0;
    }

    // Function to print the matrix
    void print() const {
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
};

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


// Wrapper function to diagonalize real symmetric matrices
inline std::tuple<Vector, DenseMatrix> solve_eigensystem(SparseMatrix A) {
    const int n = A.getNumRows();

    if (n != A.getNumCols()) {
        throw std::runtime_error("Matrix is not square");
    }

    Vector E(n);


    // over her create a dense matrix from the sparse matrix 
    DenseMatrix denseA = A;
    //then get the deep copy 
    
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, denseA.access_raw_pointer(), n, E.access_raw_pointer());

    if (info != 0) {
        throw std::runtime_error("LAPACK exited with code " + std::to_string(info));
    }

    return {E, denseA};
}

int main() {
    // Create a symmetric matrix (3x3 for simplicity)
    // populate that matrix w/ values to test it 

    // Create a sparse matrix with 3 rows and 3 columns
    SparseMatrix A(3, 3);

    // Insert some non-zero elements
    A.insert(0, 0, 1.0); // Element at (0, 0)
    A.insert(1, 1, 2.0); // Element at (1, 1)
    A.insert(2, 2, 3.0); // Element at (2, 2)

    // Compute eigenvalues and eigenvectors
    try {
        auto [E, Q] = solve_eigensystem(A);
        std::cout << "Eigenvalues: ";
        for (int i = 0; i < E.size(); i++) {
            std::cout << E.access_raw_pointer()[i] << " ";
        }
        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
