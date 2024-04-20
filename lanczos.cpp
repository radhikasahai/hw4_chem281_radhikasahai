#include <iostream>
#include <armadillo>

void lanczos(const arma::sp_mat& H, int k) {
    int n = H.n_rows;

    // Lanczos vectors and tridiagonal matrix
    arma::mat T(k, k, arma::fill::zeros);
    arma::mat V(n, k);

    // Initialization
    arma::vec x = arma::randu<arma::vec>(n);
    arma::vec q = x / arma::norm(x);
    arma::vec r = H * q;
    double alpha = arma::dot(q, r);
    r -= alpha * q;

    // Lanczos iterations
    for (int j = 0; j < k; j++) {
        V.col(j) = q;

        if (j > 0)
            r -= arma::dot(V.col(j - 1), r) * V.col(j - 1);

        double beta = arma::norm(r);
        T(j, j) = alpha;

        if (j < k - 1) {
            T(j + 1, j) = T(j, j + 1) = beta;
            q = r / beta;
            r = H * q - beta * V.col(j);
            alpha = arma::dot(q, r);
        }
    }

    // Diagonalize tridiagonal matrix
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, T);

    std::cout << "Tridiagonal matrix T: " << std::endl << T << std::endl;
    std::cout << "Eigenvalues of T: " << std::endl << eigval << std::endl;
}

int main() {
    // Example sparse matrix
    arma::sp_mat H = {{4, -1, 0}, {-1, 4, -1}, {0, -1, 4}};

    // Sparsity fraction
    double sparsity = 0.3; // Example sparsity fraction

    // Generate a random sparse matrix with specified sparsity
    H.randn(H.n_rows, H.n_cols);
    H.elem(arma::find(arma::randu<arma::sp_mat>(H.n_rows, H.n_cols) > sparsity)).zeros();

    int k = 3; // Subspace size
    lanczos(H, k);

    return 0;
}
