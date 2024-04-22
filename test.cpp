#include <iostream>
#include <vector>
#include <chrono> // For timing
#include "solvers.h"
#include "utils.h"

void run_test(int matrix_size, double sparsity) {
    // Create a sparse matrix
    SparseMatrix H(matrix_size, matrix_size);
    H.fillRandom(); // Fill with random values
    int k = 3; // Subspace size

    // Perform Davidson method
    std::vector<double> davidson_eigenvalues;
    auto start_time_davidson = std::chrono::high_resolution_clock::now();
    davidson(H, k, 50, davidson_eigenvalues);
    auto end_time_davidson = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_davidson = end_time_davidson - start_time_davidson;

    // Perform Lanczos method
    std::vector<double> lanczos_eigenvalues;
    auto start_time_lanczos = std::chrono::high_resolution_clock::now();
    lanczos(H, k, lanczos_eigenvalues);
    auto end_time_lanczos = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_lanczos = end_time_lanczos - start_time_lanczos;

    // Output results
    std::cout << "Matrix Size: " << matrix_size << ", Sparsity: " << sparsity << std::endl;
    std::cout << "Davidson: Time = " << duration_davidson.count() << " seconds, Iterations = " << davidson_eigenvalues.size() << std::endl;
    std::cout << "Lanczos: Time = " << duration_lanczos.count() << " seconds, Iterations = " << lanczos_eigenvalues.size() << std::endl;
}

int main() {
    // Define test parameters
    std::vector<int> matrix_sizes = {100, 200};
    std::vector<double> sparsity_levels = {0.1, 0.3};

    // Run tests for each combination of parameters
    for (int size : matrix_sizes) {
        for (double sparsity : sparsity_levels) {
            run_test(size, sparsity);
            std::cout << std::endl;
        }
    }

    return 0;
}
