#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>

class Vector {
private:
    std::vector<double> data;

public:
    // Constructor
    Vector();

    // Constructor with initial size and value
    Vector(int size, double value = 0.0);

    // Destructor
    ~Vector();

    // Function to get the size of the vector
    int size() const;

    // Function to access data by index
    double& operator[](int index);

    // Function to access data by index (const version)
    const double& operator[](int index) const;

    // Function to resize the vector
    void resize(int newSize, double value = 0.0);

    // Function to clear the vector
    void clear();

    // Function to compute the dot product with another vector
    double dot(const Vector& other) const;

    // Function to print the vector
    void print() const;
};

#endif // VECTOR_H
