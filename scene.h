#pragma once

#include <material.h>
#include <vector.h>
#include <object.h>
#include <light.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
#include <optional>

#define GETINDEX(vector, id) (id < 0 ? vector.size() + id : id - 1)

class Scene {
public:
    std::vector<Object> objects_;
    std::vector<SphereObject> sphere_objects_;
    std::vector<Light> lights_;
    std::unordered_map<std::string, Material> materials_;
    Scene() = default;
    Scene(std::vector<Object> objects, std::vector<SphereObject> sphere_objects,
          std::vector<Light> lights, std::unordered_map<std::string, Material> materials)
        : objects_(objects),
          sphere_objects_(sphere_objects),
          lights_(lights),
          materials_(materials) {
    }
    const std::vector<Object>& GetObjects() const {
        return objects_;
    }
    const std::vector<SphereObject>& GetSphereObjects() const {
        return sphere_objects_;
    }
    const std::vector<Light>& GetLights() const {
        return lights_;
    }
    void SetObjects(std::vector<Object>& objects) {
        objects_ = objects;
    }
    const std::unordered_map<std::string, Material>& GetMaterials() const {
        return materials_;
    }
    void AddObject(const Object& object) {
        objects_.push_back(object);
    }
    void AddSphereObject(const SphereObject& sphere_object) {
        sphere_objects_.push_back(sphere_object);
    }
    void AddLight(const Light& light) {
        lights_.push_back(light);
    }
    void AddMaterial(const std::string& name, const Material& material) {
        materials_.insert({name, material});
    }
};

std::unordered_map<std::string, Material> ReadMaterials(const std::filesystem::path& path) {
    std::unordered_map<std::string, Material> result;
    std::ifstream file;
    file.open(path);
    std::string line, name;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string command;
        ss >> command;
        if (command.empty()) {
            continue;
        }
        if (command == "newmtl") {
            ss >> name;
            Material material;
            material.name = name;
            material.albedo = {1, 0, 0};
            result[name] = material;
        } else if (!name.empty()) {
            if (command == "Ka") {
                double r, g, b;
                ss >> r >> g >> b;
                result[name].ambient_color = Vector(r, g, b);
            } else if (command == "Kd") {
                double r, g, b;
                ss >> r >> g >> b;
                result[name].diffuse_color = Vector(r, g, b);
            } else if (command == "Ks") {
                double r, g, b;
                ss >> r >> g >> b;
                result[name].specular_color = Vector(r, g, b);
            } else if (command == "Ke") {
                double r, g, b;
                ss >> r >> g >> b;
                result[name].intensity = Vector(r, g, b);
            } else if (command == "Ns") {
                double ns;
                ss >> ns;
                result[name].specular_exponent = ns;
            } else if (command == "Ni") {
                double ni;
                ss >> ni;
                result[name].refraction_index = ni;
            } else if (command == "al") {
                double cx, cy, cz;
                ss >> cx >> cy >> cz;
                result[name].albedo = {cx, cy, cz};
            }
        }
    }
    file.close();
    return result;
}

std::pair<int, int> ParseString(const std::string& s) {
    std::pair<int, int> result;
    std::istringstream ss(s);
    std::string token;

    // Split the string by "//" delimiter
    while (std::getline(ss, token, '/')) {
        if (!token.empty()) {
            if (result.first == 0) {
                result.first = std::stoi(token);
            } else {
                result.second = std::stoi(token);
            }
        }
    }

    return result;
}
auto split = [](const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream token_stream(s);
    while (std::getline(token_stream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
};

std::string RemoveSpaces(std::string& s) {
    int last = s.size() - 1;
    while (last >= 0 && s[last] == ' ') {
        --last;
    }
    return s.substr(0, last + 1);
}

void Parse(std::stringstream& ss, Material* current_material, const std::vector<Vector> vertices,
           const std::vector<Vector> normals, std::vector<Object>& objs) {
    if (ss.str().empty()) {
        return;
    }
    std::string line = ss.str();
    line.erase(0, 2);
    line = RemoveSpaces(line);
    std::string v;
    std::vector<std::pair<Vector, std::optional<Vector>>> face_v;
    std::vector<std::pair<int, int>> pairs;
    std::vector<std::string> tokens = split(line, ' ');
    std::vector<std::pair<int, int>> index;
    for (auto elem : tokens) {
        auto indices = ParseString(elem);
        std::pair<Vector, std::optional<Vector>> result;
        result.first = vertices[GETINDEX(vertices, indices.first)];
        if (GETINDEX(normals, indices.second) >= 0 &&
            GETINDEX(normals, indices.second) < normals.size()) {
            result.second = normals[GETINDEX(normals, indices.second)];
        }
        face_v.push_back(result);
    }
    for (size_t i = 0; i + 2 < face_v.size(); ++i) {
        if (face_v[0].second != std::nullopt) {
            std::array<Vector, 3> norms = {face_v[0].second.value(), face_v[i + 1].second.value(),
                                           face_v[i + 2].second.value()};
            Triangle triangle = {face_v[0].first, face_v[i + 1].first, face_v[i + 2].first};
            Object tmp = {current_material, triangle, norms};
            objs.push_back(tmp);
        } else {
            Triangle triangle = {face_v[0].first, face_v[i + 1].first, face_v[i + 2].first};
            Object tmp = {current_material, triangle, std::nullopt};
            objs.push_back(tmp);
        }
    }
}
Scene ReadScene(const std::filesystem::path& path) {
    Scene result;
    std::vector<Object> objects;
    std::vector<Light> lights;
    std::vector<SphereObject> sphere_objects;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::ifstream file;
    file.open(path);
    std::string line, current_material;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::stringstream ss(line);
        std::string command;
        ss >> command;
        if (command == "v") {
            double x, y, z;
            ss >> x >> y >> z;
            vertices.push_back(Vector(x, y, z));
        } else if (command == "vn") {
            double x, y, z;
            ss >> x >> y >> z;
            normals.push_back(Vector(x, y, z));
        } else if (command == "S") {
            double x, y, z, r;
            ss >> x >> y >> z >> r;
            Sphere sphere = Sphere(Vector(x, y, z), r);
            SphereObject sphere_object = {&result.materials_[current_material], sphere};
            sphere_objects.push_back(sphere_object);
        } else if (command == "f") {
            Parse(ss, &result.materials_[current_material], vertices, normals, objects);
        } else if (command == "mtllib") {
            std::string filename;
            std::string path_s = static_cast<std::string>(path);
            size_t pos = path_s.find_last_of('/');
            ss >> filename;
            result.materials_ =
                ReadMaterials(static_cast<std::string>(path_s.substr(0, pos + 1)) + filename);
        } else if (command == "usemtl") {
            ss >> current_material;
        } else if (command == "P") {
            double x, y, z, r, g, b;
            ss >> x >> y >> z >> r >> g >> b;
            lights.push_back(Light(Vector(x, y, z), Vector(r, g, b)));
        } else if (command == "illum") {
        } else {
            continue;
        }
    }
    file.close();

    result.lights_ = lights;
    result.objects_ = objects;
    result.SetObjects(objects);
    result.sphere_objects_ = sphere_objects;
    return result;
}
