#include <iostream>
#include <vector>
#include <tuple>
#include <lapacke.h>  
// g++ -o lapacke lapacke.cpp -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas
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
inline std::tuple<sjc::vec, sjc::mat> solve_eigensystem(sjc::mat A) {
    const int n = A.get_rows();
    if (n != A.get_cols()) {
        throw std::runtime_error("Matrix is not square");
    }

    sjc::vec E(n);
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A.access_raw_pointer(), n, E.access_raw_pointer());

    if (info != 0) {
        throw std::runtime_error("LAPACK exited with code " + std::to_string(info));
    }

    return {E, A};
}

int main() {
    // Create a symmetric matrix (3x3 for simplicity)
    sjc::mat A(3, 3);
    std::vector<double> sym_matrix = {1, 2, 3, 2, 5, 4, 3, 4, 6};
    std::copy(sym_matrix.begin(), sym_matrix.end(), A.access_raw_pointer());

    // Compute eigenvalues and eigenvectors
    try {
        auto [E, Q] = solve_eigensystem(A);
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
