
#include <valarray>
#include <iostream>
#include <initializer_list> // For std::initializer_list

template <typename T>
class my_vector : public std::valarray<T> {
public:
    // Constructor from size and initial value (matches one of the valarray constructors)
    my_vector(size_t n, T val = T()) : std::valarray<T>(val, n) {}

    // Constructor from an initializer list
    my_vector(std::initializer_list<T> il) : std::valarray<T>(il.begin(), il.size()) {}

    // Copy constructor
    my_vector(const std::valarray<T>& other) : std::valarray<T>(other) {}

    // Example of an additional method that could be useful
    T sum() const {
        T sum = 0;
        for (size_t i = 0; i < this->size(); ++i) {
            sum += (*this)[i];
        }
        return sum;
    }

    // Example method to display the contents of the vector
    void display() const {
        for (size_t i = 0; i < this->size(); ++i) {
            std::cout << (*this)[i] << ' ';
        }
        std::cout << std::endl;
    }
};

int main() {
    // Example usage with an initializer list
    my_vector<double> vec = {1.0, 2.0, 3.0, 4.0};
    vec.display();
    std::cout << "Sum: " << vec.sum() << std::endl;

    return 0;
}
