#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>
#include <iostream>

class Vector {
public:
    // initialization

    Vector() : data_(std::array<double, 3>{0, 0, 0}) {
    }
    Vector(double x, double y, double z) : data_(std::array<double, 3>{x, y, z}) {
    }
    // access

    double& operator[](size_t ind) {
        return data_[ind];
    }
    double operator[](size_t ind) const {
        return data_[ind];
    }
    // vector operations
    Vector& operator+=(const Vector& other) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] += other[i];
        }
        return *this;
    }
    Vector operator+(const Vector& other) {
        Vector temp = *this;
        temp += other;
        return temp;
    }
    const Vector operator+(const Vector& other) const {
        const Vector temp = {data_[0] + other[0], data_[1] + other[1], data_[2] + other[2]};
        return temp;
    }
    Vector& operator-=(const Vector& other) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] -= other[i];
        }
        return *this;
    }
    Vector operator-(const Vector& other) const {
        Vector temp = *this;
        temp -= other;
        return temp;
    }
    // scalar operations
    template <typename Numeric>
    Vector& operator+=(Numeric num) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] += num;
        }
        return *this;
    }
    template <typename Numeric>
    Vector operator+(Numeric num) const {
        Vector temp = *this;
        temp += num;
        return temp;
    }
    template <typename Numeric>
    friend Vector operator+(Numeric num, const Vector& self) {
        return Vector(self[0] + num, self[1] + num, self[2] + num);
    }
    template <typename Numeric>
    Vector& operator-=(Numeric num) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] -= num;
        }
        return *this;
    }
    template <typename Numeric>
    Vector operator-(Numeric num) const {
        Vector temp = *this;
        temp -= num;
        return temp;
    }
    template <typename Numeric>
    friend Vector operator-(Numeric num, const Vector& self) {
        return Vector(self[0] - num, self[1] - num, self[2] - num);
    }
    template <typename Numeric>
    Vector& operator*=(Numeric num) {
        for (size_t i = 0; i < 3; ++i) {
            data_[i] *= num;
        }
        return *this;
    }
    template <typename Numeric>
    Vector operator*(Numeric num) {
        Vector temp = *this;
        temp *= num;
        return temp;
    }
    template <typename Numeric>
    const Vector operator*(Numeric num) const {
        const Vector temp = {data_[0] * num, data_[1] * num, data_[2] * num};
        return temp;
    }
    template <typename Numeric>
    friend Vector operator*(Numeric num, const Vector& self) {
        return Vector(self[0] * num, self[1] * num, self[2] * num);
    }
    // Normalization
    double Norm() const {
        return std::sqrt(data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2]);
    }
    void Normalize() {
        double norm = this->Norm();
        for (size_t i = 0; i < 3; ++i) {
            data_[i] = data_[i] / norm;
        }
    }
    bool EqualToZero() {
        if (std::sqrt(data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2]) < 1e-9) {
            return true;
        } else {
            return false;
        }
    }

private:
    std::array<double, 3> data_;
};

inline double DotProduct(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vector CrossProduct(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
inline double Length(const Vector& v) {
    return v.Norm();
}

inline Vector HadamardProduct(const Vector& lhs, const Vector& rhs) {
    return {lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2]};
}
inline Vector HadamardDivision(const Vector& lhs, const Vector& rhs) {
    return {lhs[0] / rhs[0], lhs[1] / rhs[1], lhs[2] / rhs[2]};
}
