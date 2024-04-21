#include <iostream>
#include <vector>
#include <tuple>
#include <lapacke.h>  
#include "sparseMatrix.h" 
// g++ -o lapacke1 lapacke1.cpp -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas

namespace sjc {
    // Simple vector class for eigenvalues
    class vec {
        std::vector<double> data;
    public:
        vec(size_t n) : data(n) {}
        double* access_raw_pointer() { return data.data(); }
        size_t size() const { return data.size(); }
    };

    // Matrix class for dense matrices
    class mat {
        std::vector<double> data;
        size_t rows, cols;
    public:
        mat(size_t rows, size_t cols) : data(rows * cols), rows(rows), cols(cols) {}
        double* access_raw_pointer() { return data.data(); }
        size_t get_rows() const { return rows; }
        size_t get_cols() const { return cols; }
    };
}

// Wrapper function to diagonalize real symmetric matrices
inline std::tuple<sjc::vec, sjc::mat> solve_eigensystem(SparseMatrix& A) {
    const int n = A.getNumRows();
    if (n != A.getNumCols()) {
        throw std::runtime_error("Matrix is not square");
    }

    // Create a dense matrix from the sparse matrix
    sjc::mat dense_A(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dense_A.access_raw_pointer()[i * n + j] = A(i, j);
        }
    }

    sjc::vec E(n);
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, dense_A.access_raw_pointer(), n, E.access_raw_pointer());

    if (info != 0) {
        throw std::runtime_error("LAPACK exited with code " + std::to_string(info));
    }

    return {E, dense_A};
}

int main() {
    // Create a sparse matrix with 3 rows and 3 columns
    SparseMatrix sparse(3, 3);

    // Insert some non-zero elements
    sparse.insert(0, 0, 1.0); // Element at (0, 0)
    sparse.insert(1, 1, 2.0); // Element at (1, 1)
    sparse.insert(2, 2, 3.0); // Element at (2, 2)

    // Compute eigenvalues and eigenvectors
    try {
        auto [E, Q] = solve_eigensystem(sparse);
        std::cout << "Eigenvalues: ";
        for (size_t i = 0; i < E.size(); i++) {
            std::cout << E.access_raw_pointer()[i] << " ";
        }
        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    return 0;
}
