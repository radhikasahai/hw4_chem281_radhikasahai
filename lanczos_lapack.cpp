#include <iostream>
#include <vector>
#include <tuple>
#include <lapacke.h>  
// g++ -o lanczos_lapack lanczos_lapack.cpp -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas

class Vector {
private:
    std::vector<double> data;

public:
    // Constructor
    Vector() {}

    // Constructor with initial size and value
    Vector(int size, double value = 0.0) : data(size, value) {}

    // Destructor
    ~Vector() {}

    double* access_raw_pointer() { return data.data(); }

    // Function to get the size of the vector
    int size() const {
        return data.size();
    }

    // Function to access data by index
    double& operator[](int index) {
        return data[index];
    }

    // Function to access data by index (const version)
    const double& operator[](int index) const {
        return data[index];
    }

    // Function to resize the vector
    void resize(int newSize, double value = 0.0) {
        data.resize(newSize, value);
    }

    // Function to clear the vector
    void clear() {
        data.clear();
    }

    void randu() {
        for (int i = 0; i < data.size(); i++) {
            double randu = static_cast<double>(rand()) / RAND_MAX;
            data[i] = randu;
        }
    }

    void multiply(double val) {
        for (int i = 0; i < data.size(); i++) {
            data[i] *= val;
        }
    }


    void divide(double val) {
        for (int i = 0; i < data.size(); i++) {
            data[i] /= val;
        }
    }
    
    void add(double val) {
        for (int i = 0; i < data.size(); i++) {
            data[i] += val;
        }
    }
    
    void subtract(double val) {
        for (int i = 0; i < data.size(); i++) {
            data[i] -= val;
        }
    }
    
    void subtractVectors(const Vector& other) {
        for (int i = 0; i < data.size(); ++i) {
            data[i] -= other[i];
        }
    }

    void norm() {
        double sum = 0.0;
        for (int i = 0; i < data.size(); i++) {
            sum += data[i];
        }
        double magnitude = std::sqrt(sum);
        this->divide(magnitude);
    }

    // Function to compute the dot product with another vector
    double dot(const Vector& other) const {
        double result = 0.0;
        for (int i = 0; i < data.size(); ++i) {
            result += data[i] * other[i];
        }
        return result;
    }

    // Function to print the vector
    void print() const {
        for (int i = 0; i < data.size(); ++i) {
            std::cout << data[i] << " ";
        }
        std::cout << std::endl;
    }
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
    Vector apply(const Vector vector) const { // vector of dot product of input against each row 
        //we have to assume that "vector" is of length of num rows or cols

        Vector result(numRows);
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

void lanczos(SparseMatrix &H, int k) { //reference to H matrix
    int n  = H.getNumRows();

    // Lanczos vectors and tridiagonal matrix
    SparseMatrix T(k,k);
    SparseMatrix V(n, k);

    // Initialization
    Vector x(n); 
    x.randu(); //x is a list of n random nums in a vec 
    x.norm(); //this is q
    
    Vector r = H.apply(x);

    double alpha = x.dot(r);
    x.multiply(alpha);

    r.subtractVectors(x);

    // //Lanczos Iterations 
    // for (int j = 0; j < k; j++) {

    //     V.col(j) = q;

    //     if (j > 0)
    //         r -= arma::dot(V.col(j - 1), r) * V.col(j - 1);

    //     double beta = arma::norm(r);
    //     T(j, j) = alpha;

    //     if (j < k - 1) {
    //         T(j + 1, j) = T(j, j + 1) = beta;
    //         q = r / beta;
    //         r = H * q - beta * V.col(j);
    //         alpha = arma::dot(q, r);
    //     }
    // }
}   

int main() {

    // Create a sparse matrix with 3 rows and 3 columns
    SparseMatrix H(3, 3);

    // Insert some non-zero elements
    H.insert(0, 0, 1.0); // Element at (0, 0)
    H.insert(1, 1, 2.0); // Element at (1, 1)
    H.insert(2, 2, 3.0); // Element at (2, 2)

    // // Compute eigenvalues and eigenvectors
    // try {
    //     auto [E, Q] = solve_eigensystem(A);
    //     std::cout << "Eigenvalues: ";
    //     for (int i = 0; i < E.size(); i++) {
    //         std::cout << E.access_raw_pointer()[i] << " ";
    //     }
    //     std::cout << std::endl;
    // } catch (const std::exception& e) {
    //     std::cerr << "Error: " << e.what() << std::endl;
    // }

    // H.randn(H.n_rows, H.n_cols);
    // H.elem(arma::find(arma::randu<arma::sp_mat>(H.n_rows, H.n_cols) > sparsity)).zeros();

    int k = 3; // Subspace size
    lanczos(H, k);
    return 0;
}
