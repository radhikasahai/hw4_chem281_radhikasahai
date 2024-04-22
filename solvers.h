#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <lapacke.h>  
#include "utils.h"

// g++ -o solvers solvers.cpp -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas

void davidson(SparseMatrix &H, int k, int max_iter, std::vector<double> &eigenvalues) {
    int n  = H.getNumRows();

    SparseMatrix B(n,k);
    B.fillRandom(); //B w/ random values 

    for (int iter = 0; iter < max_iter; iter++) {
        // Orthogonalize
        for (int j = 0; j < k; j++) {
            Vector bCol(n);
            Vector bOrtho(n);

            //orthogonalize 

            //calc the ortho vector for a specific colomun in B against the entire B matrix 
            //normalize that vector 
            //set that vector as the new column for B at j 

            for (int i = 0; i < n; i++) { 
                //first get B.col(j)
                bCol[i] = B(i,j); 
            }

            bOrtho = B.orthogonalize(bCol); //calculate bOrtho
            bOrtho.norm(); //normalize bCol
            
            //update B w/ the new ortho
            for (int i = 0; i < n; i++) {
                // basically put bDense into B 
                B.insert(i,j,bOrtho[i]);
            }
        }

        DenseMatrix hDense = H;
        DenseMatrix bDense = B;

        DenseMatrix HB(n,k);

        HB = hDense.matrixMultiply(bDense); //make a matrix mult
        DenseMatrix bT = bDense.transpose();

        DenseMatrix K(n,k);
        K = bT.matrixMultiply(HB);

        auto [E, Q] = solve_eigensystem(K); 
        Q.absolute(); //positive definite 
        // Q.print();

        DenseMatrix bNew(n,k);

        bNew = bDense.matrixMultiply(Q); //multiply B by the eigenvectors 

        //subtract bDense from bNew
        bNew = bNew.subtractMatricies(bDense);

        double mag = bNew.matrixMagnitude(); 

        std::cout << mag << std::endl;
        // Check for convergence
        if (mag < 1e-3) {
            std::cout << "Davidson eigenvalues" << std::endl;
            E.print();
            std::cout << "Davidson eigenvectors" << std::endl;
            Q.print();
            std::cout << "Number of iterations: " << iter << std::endl;
            break;
        }
        //update B w/ the new bDense 
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < k; j++) {
                // basically put bDense into B 
                B.insert(i,j,bDense(i,j));
            }
        }
    }
}

void lanczos(SparseMatrix &H, int k, std::vector<double> &eigenvalues) { //reference to H matrix
    int n  = H.getNumRows();

    // Lanczos vectors and tridiagonal matrix
    SparseMatrix T(k,k);
    SparseMatrix V(n, k);

    // Initialization
    Vector x(n); 
    x.randu(); //x is a list of n random nums in a vec 
    x.norm(); //this is q

    Vector q = x; 

    Vector r = H.apply(x);

    double alpha = x.dot(r);
    x.multiply(alpha);

    r.subtractVectors(x);

    //Lanczos Iterations 
    for (int j = 0; j < k; j++) {
        Vector vCol(n);

        // here we need to set every value of the column j as the value in q 
        for (int i = 0; i < n; i++) { //V.col(j) = q;
            V.insert(i,j,q[i]); 
            vCol[i] = V(i,j-1); // get V.col(j-1);
        }

        if (j > 0) {
            double dot = vCol.dot(r);
            vCol.multiply(dot);
            r.subtractVectors(vCol);
        }

        double beta = r.magnitude();
        T.insert(j, j, alpha); 


        if (j < k - 1) {
            T.insert(j + 1, j, beta); 
            T.insert(j, j + 1, beta); 

            r.divide(beta);
            q = r;

            vCol.multiply(beta);
            r = H.apply(q);
            r.subtractVectors(vCol);

            alpha = q.dot(r);
        }
    }

    //Diagonilze Tridigonal Matrix thingy 

    // Compute eigenvalues and eigenvectors
    try {
        DenseMatrix denseT = T;
        auto [E, Q] = solve_eigensystem(denseT);
        std::cout << "Eigenvalues: ";
        for (int i = 0; i < E.size(); i++) {
            std::cout << E.access_raw_pointer()[i] << " ";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

}   

