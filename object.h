#pragma once

#include <triangle.h>
#include <material.h>
#include <sphere.h>
#include <vector.h>
#include <optional>

struct Object {
    Material* material = nullptr;
    Triangle polygon;
    std::optional<std::array<Vector, 3>> normals_;
    const Vector* GetNormal(size_t index) const {
        if (normals_.has_value()) {
            return &normals_.value()[index];
        } else {
            return nullptr;
        }
    }
    bool HasNormal() const {
        return normals_.has_value();
    }
    void SetNormal(Vector v, size_t ind) {
        normals_.value()[ind] = v;
    }
};

struct SphereObject {
    const Material* material = nullptr;
    Sphere sphere;
};
