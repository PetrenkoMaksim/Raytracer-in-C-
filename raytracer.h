#pragma once

#include <cmath>
#include <filesystem>
#include <geometry.h>
#include <image.h>
#include <intersection.h>
#include <light.h>
#include <material.h>
#include <matrix.h>
#include <object.h>
#include <optional>
#include <options/camera_options.h>
#include <options/render_options.h>
#include <ray.h>
#include <scene.h>
#include <sphere.h>
#include <triangle.h>
#include <utility>
#include "options/render_options.h"
#include <vector.h>
#include <cmath>
#include <vector>

const Vector kK = {0, 0, 1};
const Vector kI = {1, 0, 0};
const Vector kJ = {0, 1, 0};

std::pair<std::optional<Intersection>, std::optional<Material>> GetNearestIntersection(
    const Scene &scene, const Ray ray) {
    std::pair<std::optional<Intersection>, std::optional<Material>> illumination = {std::nullopt,
                                                                                    std::nullopt};
    double min_distance = -1;
    bool not_seen = true;
    for (Object object : scene.GetObjects()) {
        auto intersection = GetIntersection(ray, object.polygon);
        if (intersection) {
            if (not_seen || intersection->GetDistance() < min_distance) {
                not_seen = false;
                min_distance = intersection->GetDistance();
                illumination.second = *object.material;
                illumination.first = intersection;
                if (object.HasNormal()) {
                    Vector coord =
                        GetBarycentricCoords(object.polygon, illumination.first->GetPosition());
                    Vector normal(coord[0] * (*object.GetNormal(0)) +
                                  coord[1] * (*object.GetNormal(1)) +
                                  coord[2] * (*object.GetNormal(2)));
                    normal.Normalize();
                    if (DotProduct(ray.GetDirection(), normal) > 0) {
                        illumination.first->SetNormal(-1 * normal);
                    } else {
                        illumination.first->SetNormal(normal);
                    }
                }
            }
        }
    }
    for (SphereObject object : scene.GetSphereObjects()) {
        auto intersection = GetIntersection(ray, object.sphere);
        if (intersection) {
            if (not_seen || intersection->GetDistance() < min_distance) {
                not_seen = false;
                min_distance = intersection->GetDistance();
                illumination.second = *object.material;
                illumination.first = intersection;
            }
        }
    }
    if (not_seen || min_distance == -1) {
        return {};
    }
    return illumination;
}

bool Zero(Vector v) {
    if (std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) < 1e-9) {
        return true;
    } else {
        return false;
    }
}

Vector GetV(const Vector &fv) {
    if (!Zero(CrossProduct(kJ, fv))) {
        return kJ;
    } else if (DotProduct(kJ, fv) < 0) {
        return kK;
    } else {
        return -1. * kK;
    }
}

Ray GetRay(const CameraOptions &camera_options, int x, int y) {
    double xx = (2 * (x + 0.5) / static_cast<double>(camera_options.screen_width) - 1) *
                tan(camera_options.fov / 2.) * camera_options.screen_width /
                static_cast<double>(camera_options.screen_height),
           yy = -(2 * (y + 0.5) / static_cast<double>(camera_options.screen_height) - 1) *
                tan(camera_options.fov / 2.),
           zz = -1;

    Vector dir = {xx, yy, zz};
    dir.Normalize();
    Vector fv = Vector(camera_options.look_from) - Vector(camera_options.look_to);
    fv.Normalize();

    Vector tmp = GetV(fv);
    tmp.Normalize();

    Vector right = CrossProduct(tmp, fv);
    right.Normalize();
    Vector up = CrossProduct(fv, right);
    up.Normalize();

    Mat transform(right, up, fv, camera_options.look_from);

    return Ray(MatrixDotProduct({0, 0, 0}, transform),
               MatrixDotProduct(dir, transform) - MatrixDotProduct({0, 0, 0}, transform));
}

enum class RefractOptions { Refract, NotRefract };

std::pair<Vector, bool> GetColor(const Scene &scene, const Ray &ray,
                                 const RenderOptions &render_options, double mx_distance,
                                 int recursion_depth = 0,
                                 RefractOptions refraction_type = RefractOptions::NotRefract);

Vector GetIlluminationBase(
    const Ray &ray, const std::pair<std::optional<Intersection>, std::optional<Material>> &nearest,
    const Scene &scene) {
    Material material = nearest.second.value();
    Intersection intersection = nearest.first.value();
    Vector illumination(0., 0., 0.);
    illumination += material.ambient_color;
    illumination += material.intensity;
    for (auto light : scene.GetLights()) {
        auto direction = light.position - intersection.GetPosition();
        Ray ray_tmp = {1e-9 * intersection.GetNormal() + intersection.GetPosition(), direction};
        auto obstacle = GetNearestIntersection(scene, ray_tmp).first;
        if (obstacle.has_value() && obstacle->GetDistance() < Length(direction)) {
            continue;
        }

        illumination +=
            material.albedo[0] *
            (HadamardProduct(material.diffuse_color, light.intensity) *
                 std::max(0., DotProduct(ray_tmp.GetDirection(), intersection.GetNormal())) +
             HadamardProduct(material.specular_color, light.intensity) *
                 pow(std::max(0., DotProduct(ray.GetDirection(),
                                             -1 * Reflect(-1 * ray_tmp.GetDirection(),
                                                          intersection.GetNormal()))),
                     material.specular_exponent));
    }

    return illumination;
}

void GetRefract(const Ray &ray,
                const std::pair<std::optional<Intersection>, std::optional<Material>> &nearest,
                const Scene &scene, int curr_depth, RefractOptions refraction_type,
                const RenderOptions &render_options, Vector &illumination) {
    double ref_index =
        ((refraction_type == RefractOptions::Refract) ? 1. / nearest.second->refraction_index
                                                      : nearest.second->refraction_index);
    if (Refract(ray.GetDirection(), nearest.first.value().GetNormal(), ref_index) &&
        std::abs(nearest.second->albedo[2]) > 1e-9) {
        Ray refract_ray =
            Ray(-1e-9 * nearest.first.value().GetNormal() + nearest.first.value().GetPosition(),
                Refract(ray.GetDirection(), nearest.first.value().GetNormal(), ref_index).value());
        if (refraction_type == RefractOptions::NotRefract) {
            illumination = illumination + GetColor(scene, refract_ray, render_options, -1.,
                                                   curr_depth + 1, RefractOptions::Refract)
                                              .first;
        } else {
            illumination = illumination + nearest.second->albedo[2] *
                                              GetColor(scene, refract_ray, render_options, -1.,
                                                       curr_depth + 1, RefractOptions::NotRefract)
                                                  .first;
        }
    }
};

void GetReflect(const Ray &ray,
                const std::pair<std::optional<Intersection>, std::optional<Material>> &nearest,
                const Scene &scene, int curr_depth, RefractOptions refraction_type,
                const RenderOptions &render_options, Vector &illumintaion) {
    if (refraction_type == RefractOptions::NotRefract &&
        std::abs(nearest.second->albedo[1]) > 1e-9) {
        Ray reflect_ray =
            Ray(1e-9 * nearest.first.value().GetNormal() + nearest.first.value().GetPosition(),
                Reflect(ray.GetDirection(), nearest.first.value().GetNormal()));
        illumintaion +=
            nearest.second->albedo[1] * GetColor(scene, reflect_ray, render_options, -1.,
                                                 curr_depth + 1, RefractOptions::NotRefract)
                                            .first;
    }
};

Vector GetIlluminationComp(
    const Ray &ray, const std::pair<std::optional<Intersection>, std::optional<Material>> &nearest,
    const Scene &scene, int curr_depth, RefractOptions refraction_type,
    const RenderOptions &render_options) {
    Vector illumination = {0, 0, 0};
    GetReflect(ray, nearest, scene, curr_depth, refraction_type, render_options, illumination);
    GetRefract(ray, nearest, scene, curr_depth, refraction_type, render_options, illumination);
    return illumination;
}

std::pair<Vector, bool> GetColor(const Scene &scene, const Ray &ray,
                                 const RenderOptions &render_options, double mx_distance,
                                 int recursion_depth, RefractOptions refraction_type) {

    Vector rgb(0., 0., 0.);
    auto nearest = GetNearestIntersection(scene, ray);
    mx_distance = std::max(mx_distance, nearest.first->GetDistance());
    if (render_options.mode == RenderMode::kDepth) {
        if (nearest.first) {
            rgb = {nearest.first->GetDistance(), nearest.first->GetDistance(),
                   nearest.first->GetDistance()};
        } else {
            return std::make_pair(rgb, false);
        }
    } else if (render_options.mode == RenderMode::kNormal) {
        if (nearest.first) {
            rgb = nearest.first->GetNormal() * (1. / 2) + (1. / 2);
        } else {
            return std::make_pair(rgb, false);
        }
    } else if (render_options.mode == RenderMode::kFull) {
        if (nearest.first && recursion_depth < render_options.depth) {
            return std::make_pair(GetIlluminationBase(ray, nearest, scene) +
                                      GetIlluminationComp(ray, nearest, scene, recursion_depth,
                                                          refraction_type, render_options),
                                  true);
        } else {
            Vector rgb = {0, 0, 0};
            return std::make_pair(rgb, true);
        }
    }
    return std::make_pair(rgb, false);
}

double GetMaxDistance(const Scene &scene, const CameraOptions &camera_options) {
    double mx_distance = -1.;

    for (int i = 0; i < camera_options.screen_height; i++) {
        for (int j = 0; j < camera_options.screen_width; j++) {
            std::optional<Intersection> intersection =
                GetNearestIntersection(scene, GetRay(camera_options, i, j)).first;

            if (intersection && mx_distance < intersection->GetDistance()) {
                mx_distance = intersection->GetDistance();
            }
        }
    }
    return mx_distance;
}
template <typename T>
T MyMax(std::vector<T> l) {
    T result = l[0];
    for (auto elem : l) {
        if (elem > result) {
            result = elem;
        }
    }
    return result;
}

double GetCnst(const Scene &scene, const CameraOptions &camera_options,
               const RenderOptions &render_options, double mx_distance) {
    double cnst = -1.;
    for (int i = 0; i < camera_options.screen_height; i++) {
        for (int j = 0; j < camera_options.screen_width; j++) {
            Vector rgb =
                GetColor(scene, GetRay(camera_options, j, i), render_options, mx_distance).first;
            std::vector rgbs = {rgb[0], rgb[1], rgb[2], cnst};
            cnst = MyMax(rgbs);
        }
    }
    return cnst;
}

Vector ToneMapping(Vector rgb, double cnst) {
    return HadamardDivision(HadamardProduct(rgb, std::pow(1 / cnst, 2) * rgb + 1), 1 + rgb);
}

RGB GetRgb(int i, int j, const RenderOptions &render_options, double mx_distance,
           std::vector<std::vector<Vector>> &rgb_matrix, double cnst) {
    RGB rgb = {255, 255, 255};
    if (render_options.mode == RenderMode::kDepth) {
        if (std::abs(mx_distance) < 1e-8) {
            return rgb;
        }
        double inv_dist = 1. / mx_distance;
        int r = static_cast<int>(255 * (rgb_matrix[j][i][0] * inv_dist));
        int g = static_cast<int>(255 * (rgb_matrix[j][i][1] * inv_dist));
        int b = static_cast<int>(255 * (rgb_matrix[j][i][2] * inv_dist));
        rgb = {r, g, b};
    } else if (render_options.mode == RenderMode::kNormal) {
        int r = static_cast<int>(255 * rgb_matrix[j][i][0]);
        int g = static_cast<int>(255 * rgb_matrix[j][i][1]);
        int b = static_cast<int>(255 * rgb_matrix[j][i][2]);
        rgb = {r, g, b};
    } else if (render_options.mode == RenderMode::kFull) {
        Vector rgb_ji(rgb_matrix[j][i][0], rgb_matrix[j][i][1], rgb_matrix[j][i][2]);
        if (Length(rgb_ji) < 1e-8) {
            rgb = {0, 0, 0};
            return rgb;
        }
        rgb_ji = ToneMapping(rgb_ji, cnst);
        int r = static_cast<int>(255 * pow(rgb_ji[0], 1 / 2.2));
        int g = static_cast<int>(255 * pow(rgb_ji[1], 1 / 2.2));
        int b = static_cast<int>(255 * pow(rgb_ji[2], 1 / 2.2));
        rgb = {r, g, b};
    }
    return rgb;
}

Image Render(const std::filesystem::path &path, const CameraOptions &camera_options,
             const RenderOptions &render_options) {
    Scene scene = ReadScene(path);
    int width = camera_options.screen_width;
    int height = camera_options.screen_height;
    Image resulting_image(width, height);
    std::vector<std::vector<Vector>> rgb_matrix(height, std::vector<Vector>(width));
    double mx_distance = -1.;

    if (render_options.mode == RenderMode::kDepth) {
        mx_distance = GetMaxDistance(scene, camera_options);
    }
    double cnst = GetCnst(scene, camera_options, render_options, mx_distance);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto ray = GetRay(camera_options, i, j);
            auto color = GetColor(scene, ray, render_options, mx_distance);
            RGB rgb;
            if (color.second) {
                rgb_matrix[j][i] = color.first;
                rgb = GetRgb(i, j, render_options, mx_distance, rgb_matrix, cnst);
            } else if (render_options.mode == RenderMode::kDepth) {
                rgb = {255, 255, 255};
            } else if (render_options.mode == RenderMode::kNormal) {
                rgb = {0, 0, 0};
            }
            resulting_image.SetPixel(rgb, j, i);
        }
    }
    return resulting_image;
}
