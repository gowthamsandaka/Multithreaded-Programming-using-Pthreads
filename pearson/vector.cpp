#include <iostream>
#include <cmath>

#include "vector.hpp"

Vector::Vector()
    : size { 0 }
    , data { nullptr }
{
    // TODO: Note Implemented Yet
}

Vector::~Vector()
{
    if (data) {
        delete[] data;
    }

    size = 0;
}

Vector::Vector(unsigned size)
    : size { size }
    , data { new double[size] }
{
    // TODO: Note Implemented Yet
}

Vector::Vector(unsigned size, double* data)
    : size { size }
    , data { data }
{
    // TODO: Note Implemented Yet
}

Vector::Vector(const Vector& other)
    : Vector { other.size }
{
    for (auto i { 0 }; i < (int)size; i++) {
        data[i] = other.data[i];
    }
}

unsigned Vector::get_size() const
{
    return size;
}

double* Vector::get_data()
{
    return data;
}

double Vector::operator[](unsigned i) const
{
    return data[i];
}

double& Vector::operator[](unsigned i)
{
    return data[i];
}

double Vector::mean() const
{
    double sum { 0 };

    for (auto i { 0 }; i < (int)size; i++) {
        sum += data[i];
    }

    return sum / static_cast<double>(size);
}

double Vector::magnitude() const
{
    auto dot_prod { dot(*this) };
    return std::sqrt(dot_prod);
}

Vector Vector::operator/(double div)
{
    auto result { *this };

    for (auto i { 0 }; i < (int)size; i++) {
        result[i] /= div;
    }

    return result;
}

Vector Vector::operator-(double sub)
{
    auto result { *this };

    for (auto i { 0 }; i < (int)size; i++) {
        result[i] -= sub;
    }

    return result;
}

double Vector::dot(Vector rhs) const
{
    double result { 0 };

    for (auto i { 0 }; i < (int)size; i++) {
        result += data[i] * rhs[i];
    }

    return result;
}
