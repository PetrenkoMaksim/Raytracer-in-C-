#pragma once

#include <vector.h>
#include <cmath>

class Triangle {
public:
    Triangle(const Vector& a, const Vector& b, const Vector& c) : a_(a), b_(b), c_(c) {
    }
    const Vector& operator[](size_t ind) const {
        if (ind == 0) {
            return a_;
        } else if (ind == 1) {
            return b_;
        } else {
            return c_;
        }
    }
    double Area() const {
        return Length(CrossProduct(b_ - a_, c_ - a_)) / 2;
    }
    Vector& GetVertex(size_t ind) {
        if (ind == 0) {
            return a_;
        } else if (ind == 1) {
            return b_;
        } else {
            return c_;
        }
    }

private:
    Vector a_;
    Vector b_;
    Vector c_;
};
