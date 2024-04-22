#include <iostream>
#include <armadillo>

class DavidsonSolver {
private:
    arma::mat H;
    int k;
    int max_iter;
    double tol;

    void orthogonalize(arma::vec& v, const arma::mat& V) {
        for (int i = 0; i < V.n_cols; i++) {
            v -= (arma::dot(V.col(i), v) / arma::dot(V.col(i), V.col(i))) * V.col(i);
        }
    }

public:
    DavidsonSolver(const arma::mat& H_, int k_, int max_iter_, double tol_)
        : H(H_), k(k_), max_iter(max_iter_), tol(tol_) {}

    void solve() {
        int n = H.n_rows;
        arma::mat B(n, k, arma::fill::randn);
        arma::mat B_new(n, k, arma::fill::zeros);

        for (int iter = 0; iter < max_iter; iter++) {
            // Orthogonalize
            for (int i = 0; i < k; i++) {
                orthogonalize(B.col(i), B);
                B.col(i) = B.col(i) / arma::norm(B.col(i));
            }
    
            // Matrix-vector product
            arma::mat HB = H * B;

            // Build the subspace matrix
            arma::mat K = B.t() * HB;

            // Diagonalize the subspace matrix
            arma::vec eigenvalues;
            arma::mat eigenvectors;
            eig_sym(eigenvalues, eigenvectors, K);

            // Select Ritz vectors (eigenvectors) with smallest eigenvalues
            arma::mat new_B = B * eigenvectors.cols(0, k - 1);

            // Check for convergence
            if (arma::norm(new_B - B) < tol) {
                eigenvalues.print("Davidson eigenvalues");
                eigenvectors.print("Davidson eigenvectors");
                std::cout << "Number of iterations: " << iter << std::endl;
                break;
            }

            // Update guess basis
            B = new_B;
        }
    }
};

int main() {
    int n = 100;
    arma::mat H(n, n, arma::fill::randn);
    H = 0.5 * (H + H.t()); // Symmetrize H

    int k = 5;             // Number of eigenvalues and eigenvectors to compute
    int max_iter = 200;
    double tol = 1e-3;

    DavidsonSolver solver(H, k, max_iter, tol);
    solver.solve();

    return 0;
}
