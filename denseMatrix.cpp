

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cassert>

template <typename T>
class DenseMatrix {
private:
    std::vector<std::vector<T>> mat;
    size_t rows, cols;

public:
    // Constructors
    DenseMatrix(size_t m, size_t n) : rows(m), cols(n), mat(m, std::vector<T>(n, T())) {}

    // Access elements using bounds-checked method
    std::vector<T>& operator[](size_t i) { return mat.at(i); }
    const std::vector<T>& operator[](size_t i) const { return mat.at(i); }

    // Bounds-checked access
    std::vector<T>& at(size_t i) {
        return mat.at(i);
    }

    const std::vector<T>& at(size_t i) const {
        return mat.at(i);
    }

    // Get dimensions
    size_t get_rows() const { return rows; }
    size_t get_columns() const { return cols; }

    void resize(size_t m, size_t n) {
        if (m < rows) {
            mat.resize(m);  // Trims excess rows
        } else if (m > rows) {
            mat.resize(m, std::vector<T>(cols, T()));  // Adds new rows with existing column size
        }
        rows = m;

        for (auto& row : mat) {
            if (n < cols) {
                row.resize(n);  // Trims excess columns
            } else if (n > cols) {
                row.resize(n, T());  // Expands each row to new column size, initializing new elements
            }
        }
        cols = n;
    }
};

int main() {
    DenseMatrix<int> mat(3, 4);
    assert(mat.get_rows() == 3);
    assert(mat.get_columns() == 4);

    mat[0][0] = 1;
    mat[0][1] = 2;
    mat[1][0] = 3;
    mat[2][3] = 4;
    assert(mat[0][0] == 1);
    assert(mat[0][1] == 2);
    assert(mat[1][0] == 3);
    assert(mat[2][3] == 4);

    // Proper out-of-bounds access check using at()
    try {
        int test = mat.at(3).at(4);  // Should throw std::out_of_range
    } catch (const std::out_of_range& e) {
        std::cout << "Out-of-bounds access caught correctly: " << e.what() << std::endl;
    }

    mat.resize(5, 5);
    assert(mat.get_rows() == 5);
    assert(mat.get_columns() == 5);
    mat[4][4] = 5;
    assert(mat[4][4] == 5);

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
