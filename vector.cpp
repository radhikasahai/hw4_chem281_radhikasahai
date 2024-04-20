#include <vector>
#include <iostream>

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

// TEST


// int main() {
//     // Create a vector with default constructor
//     Vector vec1;

//     // Print size of vec1
//     std::cout << "Size of vec1: " << vec1.size() << std::endl;

//     // Resize vec1 to size 5 with value 2
//     vec1.resize(5, 2.0);

//     // Print dataa of vec1
//     std::cout << "data of vec1: ";
//     vec1.print();

//     // Access and modify data of vec1
//     vec1[0] = 1.0;
//     vec1[3] = 3.0;

//     // Print data of vec1 after modification
//     std::cout << "data of vec1 after modification: ";
//     vec1.print();

//     // Create a vector with initial size 3 and value 4
//     Vector vec2(3, 4.0);

//     // Print data of vec2
//     std::cout << "data of vec2: ";
//     vec2.print();

//     // Compute dot product of vec1 and vec2
//     double dotProduct = vec1.dot(vec2);
//     std::cout << "Dot product of vec1 and vec2: " << dotProduct << std::endl;

//     return 0;
// }
