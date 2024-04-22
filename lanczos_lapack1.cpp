#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <lapacke.h>  
// g++ -o lanczos_lapack1 lanczos_lapack1.cpp -std=c++17 -I/opt/homebrew/Cellar/lapack/3.12.0/include -L/opt/homebrew/lib -llapacke -llapack -lblas


class Vector {
private:
    std::vector<double> data;

public:
    // Constructors
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

    void multiplyInverse(double val) {
        for (int i = 0; i < data.size(); i++) {
            data[i] = val / data[i];
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

    
    void multiplyVectors(const Vector& other) {
        for (int i = 0; i < data.size(); ++i) {
            data[i] *= other[i];
        }
    }

    double magnitude() {
        double sum = 0.0;
        for (int i = 0; i < data.size(); i++) {
            sum += data[i];
        }
       return std::sqrt(std::abs(sum));
    }

    void norm() {
        this->divide(this->magnitude());
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

class SparseMatrix {
private:
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<double>> values;
    int numRows;
    int numCols;

public:
    // Constructors
    SparseMatrix() : numRows(0), numCols(0) {}

    SparseMatrix(int rows, int cols) : numRows(rows), numCols(cols) {
        indices.resize(rows);
        values.resize(rows);
    }

    // Destructor
    ~SparseMatrix() {}

    // Function to get the number of rows
    int getNumRows() const {
        return numRows;
    }

    // Function to get the number of columns
    int getNumCols() const {
        return numCols;
    }

    // Function to insert a non-zero value at a specific position
    void insert(int row, int col, double value) {
        indices[row].push_back(col);
        values[row].push_back(value);
    }

    // Function to access elements by indices
    double operator()(int row, int col) const {
        for (int i = 0; i < indices[row].size(); ++i) {
            if (indices[row][i] == col) {
                return values[row][i];
            }
        }
        return 0.0; // Return 0 if the element is zero or not found
    }

    // Function to clear the sparse matrix
    void clear() {
        indices.clear();
        values.clear();
        numRows = 0;
        numCols = 0;
    }

    // Function to print the sparse matrix
    void print() const {
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }

    double* access_raw_pointer() {
        return nullptr; // Not implemented for sparse matrices
    }
};

std::tuple<Vector, DenseMatrix> solve_eigensystem(SparseMatrix &H) {
    // Solve eigenvalue problem for symmetric matrix H using LAPACK
    int n = H.getNumRows();
    int lda = n;
    int ldz = 1;
    int info;
    Vector eigenvalues(n);
    DenseMatrix eigenvectors(n, n);

    // Temporary storage
    std::vector<double> work(3 * n - 2);
    std::vector<int> iwork(3 * n - 2);

    // Call LAPACK function
    dsyev_("V", "U", &n, reinterpret_cast<double*>(H.access_raw_pointer()), &lda, reinterpret_cast<double*>(eigenvalues.access_raw_pointer()), work.data(), &n, reinterpret_cast<double*>(eigenvectors.access_raw_pointer()), &lda, &info);

    // Ensure 'E' is declared and assigned in the 'solve_eigensystem' function
    auto [E, Q] = solve_eigensystem(T);

    // Fix the range expression for iterating over Eigenvalues
    for (double val : std::get<Vector>(E)) {
        std::cout << val << " ";
    }
        // Fill eigenvectors matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            eigenvectors(i, j) = 0.0;
        }
        eigenvectors(i, i) = 1.0;
    }

    return {eigenvalues, eigenvectors};
}

std::tuple<Vector, DenseMatrix> lanczos(SparseMatrix &H, int k) {
    try {
        SparseMatrix T(k, k);
        SparseMatrix V(H.getNumRows(), k);

        // Lanczos method implementation...

        // Compute eigenvalues and eigenvectors
        auto [E, Q] = solve_eigensystem(T);

        // Print eigenvalues
        std::cout << "Lanczos Method Eigenvalues: ";
        for (int i = 0; i < E.size(); i++) {
            std::cout << E.access_raw_pointer()[i] << " ";
        }
        std::cout << std::endl;

        return {E, Q};
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return {};
    }
}

std::tuple<Vector, DenseMatrix> davidson(SparseMatrix &H, int k, int max_iter) {
    for (int iter = 0; iter < max_iter; iter++) {
        try {
            // Davidson method iterations...

            // Compute eigenvalues and eigenvectors
            auto [E, Q] = solve_eigensystem(H);

            std::cout << "Davidson Method Eigenvalues: ";
            for (double val : E) {
                std::cout << val << " ";
            }
            std::cout << std::endl;

            std::cout << "Davidson Method Eigenvectors:" << std::endl;
            Q.print();

            return {E, Q};
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return {};
        }
    }
}

int main() {
    // Example usage of Lanczos method
    SparseMatrix H(5, 5);
    // Initialize H matrix...

    // Call Lanczos method
    auto [E_lanczos, Q_lanczos] = lanczos(H, 3);

    // Example usage of Davidson method
    // Call Davidson method
    auto [E_davidson, Q_davidson] = davidson(H, 3, 10);

    return 0;
}
