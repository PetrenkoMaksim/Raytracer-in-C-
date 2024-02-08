#pragma once

#include <vector.h>

class Intersection {
public:
    Intersection() = default;
    Intersection(const Vector& position, const Vector& normal, double distance)
        : position_(position), normal_(normal), distance_(distance) {
        normal_.Normalize();
    }

    const Vector& GetPosition() const {
        return position_;
    }
    const Vector& GetNormal() const {
        return normal_;
    }
    double GetDistance() const {
        return distance_;
    }
    void SetNormal(Vector normal) {
        normal_ = normal;
    }

private:
    Vector position_;
    Vector normal_;
    double distance_;
};
