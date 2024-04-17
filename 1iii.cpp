
#include <iostream>
#include <cassert>
#include <vector>
#include <list>
#include <stdexcept>


template <typename T>
class MyVector {
public:
    std::vector<T> data;

    MyVector(size_t n = 0, T value = T()) : data(n, value) {}

    T& operator[](size_t i) { return data[i]; }
    const T& operator[](size_t i) const { return data[i]; }

    size_t size() const { return data.size(); }
};

template <typename T>
class SparseMatrix {
private:
    struct Element {
        int col;
        T value;
    };

    std::vector<std::list<Element>> rows;
    size_t num_rows, num_cols;

public:
    SparseMatrix(size_t m, size_t n) : num_rows(m), num_cols(n), rows(m) {}

    void insert(size_t i, size_t j, T val) {
        if (i >= num_rows || j >= num_cols) throw std::out_of_range("Index out of bounds.");
        rows[i].push_back({static_cast<int>(j), val});  // Corrected line
    }

    size_t get_rows() const { return num_rows; }
    size_t get_columns() const { return num_cols; }

    MyVector<T> operator*(const MyVector<T>& vec) const {
        MyVector<T> result(num_rows, T());
        for (size_t i = 0; i < num_rows; ++i) {
            T sum = T();
            for (const auto& e : rows[i]) {
                sum += e.value * vec[e.col];
            }
            result[i] = sum;
        }
        return result;
    }
};

int main() {
    SparseMatrix<int> mat(3, 3);
    mat.insert(0, 0, 1);
    mat.insert(0, 1, 2);
    mat.insert(1, 0, 3);
    mat.insert(2, 2, 4);

    assert(mat.get_rows() == 3);
    assert(mat.get_columns() == 3);

    MyVector<int> vec(3);
    vec[0] = 1;
    vec[1] = 2;
    vec[2] = 3;

    MyVector<int> result = mat * vec;
    assert(result[0] == 5);   // 1*1 + 2*2
    assert(result[1] == 3);   // 3*1
    assert(result[2] == 12);  // 4*3

    std::cout << "All tests passed!" << std::endl;
    return 0;
}