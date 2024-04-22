#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
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

        // Matrix multiplication function
    void absolute() {
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                data[i][j] = abs(data[i][j]);
                // result(i,j) = data[i][j] - mat(i, j);
            }
        }
    }

    // Function to calculate the magnitude of the matrix
    double matrixMagnitude() const {
        double sum_of_squares = 0.0;
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                sum_of_squares += data[i][j] * data[i][j];
            }
        }
        return std::sqrt(sum_of_squares);
    }

    // Matrix multiplication function
    DenseMatrix subtractMatricies(DenseMatrix& mat) {
        DenseMatrix result(numRows, numCols);
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                result(i,j) = data[i][j] - mat(i, j);
            }
        }
        return result;
    }


        // Matrix multiplication function
    DenseMatrix matrixMultiply(const DenseMatrix& other) const {
        if (numCols != other.getNumRows()) {
            throw std::invalid_argument("Matrices cannot be multiplied: Invalid dimensions.");
        }

        DenseMatrix result(numRows, other.getNumCols());
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < other.getNumCols(); ++j) {
                for (int k = 0; k < numCols; ++k) {
                    result(i, j) += data[i][k] * other(k, j);
                }
            }
        }

        return result;
    }

    // Transpose function
    DenseMatrix transpose() const {
        DenseMatrix result(numCols, numRows);
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                result(j, i) = data[i][j];
            }
        }
        return result;
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
    void fillRandom() {
        // Create a random number generator engine
        std::random_device rd;
        std::mt19937 gen(rd()); // Mersenne Twister engine
        std::normal_distribution<double> distribution(0.0, 1.0); // Standard normal distribution

        // Fill the matrix with random numbers
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                double randn = distribution(gen); // Generate random number from standard normal distribution
                // if (randn < 0) {
                //     randn *= -1;
                // }
                this->insert(i,j,randn);
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

    // // Function to access data in the matrix
    double get(int i, int j) const {
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

    // Function to apply the sparse matrix to a vector
    Vector orthogonalize(const Vector vector) const { // vector of dot product of input against each row 
        //we have to assume that "vector" is of length of num rows or cols
        Vector result(numRows);
        result = vector; //set the result to initially be the input vector

        for (int j = 0; j < numCols; j++) {
            Vector colVec(numRows);
            for (int i = 0; i < numRows; i++) {
                colVec[i] = this->get(i,j); 

            }
            //result -= dot(result, colVec) / (dot(colVec,colVec) * colVec)
            double term1 = result.dot(colVec); //dot(result, colVec)
            double term2 = colVec.dot(colVec); //dot(colVec,colVec)

            colVec.multiply(term2); // (dot(colVec,colVec) * colVec)
            colVec.multiplyInverse(term1); //dot(result, colVec) / (dot(colVec,colVec) * colVec)

            result.subtractVectors(colVec); //result -= dot(result, colVec) / (dot(colVec,colVec) * colVec)

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
inline std::tuple<Vector, DenseMatrix> solve_eigensystem(DenseMatrix &A) {
    const int n = A.getNumRows();

    if (n != A.getNumCols()) {
        throw std::runtime_error("Matrix is not square");
    }

    Vector E(n);
    // over her create a dense matrix from the sparse matrix 
    //then get the deep copy 
    
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A.access_raw_pointer(), n, E.access_raw_pointer());

    if (info != 0) {
        throw std::runtime_error("LAPACK exited with code " + std::to_string(info));
    }
    return {E, A};
}


#endif // UTILS_H