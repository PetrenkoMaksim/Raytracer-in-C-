#pragma once

#include <vector.h>
#include <sphere.h>
#include <intersection.h>
#include <triangle.h>
#include <ray.h>

#include <optional>
#include "vector.h"

const double kEpsilon = 1e-9;

std::optional<Intersection> GetIntersection(const Ray& ray, const Sphere& sphere) {
    Vector l = sphere.GetCenter() - ray.GetOrigin();
    double tca = DotProduct(l, ray.GetDirection());
    double d2 = DotProduct(l, l) - tca * tca;
    if (d2 > sphere.GetRadius() * sphere.GetRadius()) {  // ray misses the sphere
        return {};
    }
    double thc = std::sqrt(sphere.GetRadius() * sphere.GetRadius() - d2);
    double t0 = tca - thc;
    double t1 = tca + thc;
    if (t0 < 0) {
        t0 = t1;
    }
    if (t0 < 0) {
        return {};
    }
    const Vector p1 = ray.GetOrigin() + t0 * ray.GetDirection();
    const Vector p2 = ray.GetOrigin() + t1 * ray.GetDirection();

    if (Length(ray.GetOrigin() - sphere.GetCenter()) > sphere.GetRadius()) {
        Vector normal = p1 - sphere.GetCenter();
        normal.Normalize();
        double dist = Length(ray.GetOrigin() - p1);
        return Intersection(p1, normal, dist);
    } else {
        Vector normal = (p2 - sphere.GetCenter()) * (-1);
        normal.Normalize();
        double dist = Length(ray.GetOrigin() - p2);
        return Intersection(p2, normal, dist);
    }
}

std::optional<Intersection> GetIntersection(const Ray& ray, const Triangle& triangle) {
    Vector v1, v2, temp_01, temp_2, temp_3;
    double det_1, det_2, det_3, det_4;
    v1 = triangle[1] - triangle[0];
    v2 = triangle[2] - triangle[0];
    temp_01 = CrossProduct(ray.GetDirection(), v2);
    det_1 = DotProduct(v1, temp_01);
    if (det_1 > -kEpsilon && det_1 < kEpsilon) {
        return std::nullopt;
    }
    det_2 = 1. / det_1;
    temp_2 = ray.GetOrigin() - triangle[0];
    det_3 = det_2 * DotProduct(temp_2, temp_01);
    if (det_3 < 0.0 || det_3 > 1.0) {
        return std::nullopt;
    }
    temp_3 = CrossProduct(temp_2, v1);
    det_4 = det_2 * DotProduct(ray.GetDirection(), temp_3);
    if (det_4 < 0.0 || det_3 + det_4 > 1.0) {
        return std::nullopt;
    }
    double t = det_2 * DotProduct(v2, temp_3);
    if (t > kEpsilon) {
        Vector out_intersection_point = ray.GetOrigin() + ray.GetDirection() * t;
        double dist = Length(out_intersection_point - ray.GetOrigin());
        Vector normal = CrossProduct(v1, v2);
        normal.Normalize();
        if (DotProduct(ray.GetDirection(), normal) > 0) {
            normal = -1 * normal;
        }
        return Intersection(out_intersection_point, normal, dist);
    } else {
        return {};
    }
}

Vector Reflect(const Vector& ray, const Vector& normal) {
    return ray - normal * DotProduct(normal, ray) * 2.;
}

std::optional<Vector> Refract(const Vector& ray, const Vector& normal, double eta) {
    double cosi = std::max(-1., std::min(1., DotProduct(normal, ray)));

    double etai = 1, etat = eta;
    Vector n = normal;
    if (cosi < 0) {
        cosi = -cosi;
    } else {
        std::swap(etai, etat);
        n = -1 * n;
    }
    double etaa = etai / etat;
    double k = 1 - (etaa * etaa) * (1 - (cosi * cosi));
    if (k <= 0) {
        return {};
    } else {
        return ray * etaa + n * (etaa * cosi - std::sqrt(k));
    }
}

Vector GetBarycentricCoords(const Triangle& triangle, const Vector& point) {
    return {Triangle({point, triangle[1], triangle[2]}).Area() / triangle.Area(),
            Triangle({point, triangle[0], triangle[2]}).Area() / triangle.Area(),
            Triangle({point, triangle[0], triangle[1]}).Area() / triangle.Area()};
}
