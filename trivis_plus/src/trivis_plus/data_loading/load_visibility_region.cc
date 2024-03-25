/**
 * File:   load_visibility_region.cc
 *
 * Date:   27.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/data_loading/load_visibility_region.h"

#include <fstream>
#include <iostream>

using namespace trivis_plus;
using namespace trivis_plus::data_loading;

using namespace trivis;
using namespace trivis::geom;

bool data_loading::SaveVisibilityRegion(
    const RadialVisibilityRegion &region,
    const std::string &file,
    std::stringstream *info
) {
    std::ofstream fs(file);
    if (!fs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for writing.\n";
        return false;
    }
    // todo: implement this!
    if (info) *info << "SaveVisibilityRegion is not implemented!\n";
    return false;
}

std::optional<trivis::RadialVisibilityRegion> data_loading::LoadVisibilityRegion(
    const std::string &file,
    std::stringstream *info
) {
    std::ifstream fs(file);
    if (!fs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for reading.\n";
        return std::nullopt;
    }
    // todo: implement this!
    if (info) *info << "LoadVisibilityRegion is not implemented!\n";
    return std::nullopt;
}

bool data_loading::LoadVisibilityRegion(
    const std::string &file,
    RadialVisibilityRegion &visibility_region
) noexcept(false) {
    std::ifstream ifs(file.c_str());
    if (ifs.fail()) {
        return false;
    }
    std::string token;
    while (true) {
        ifs >> token;
        if (ifs.eof()) {
            break;
        }
        if (token == "[SEED]") {
            ifs >> visibility_region.seed.x;
            ifs >> visibility_region.seed.y;
            int seed_id;
            ifs >> seed_id;
            if (seed_id != -1) {
                visibility_region.seed_id = seed_id;
            }
            continue;
        }
        if (token == "[RADIUS]") {
            double radius;
            ifs >> radius;
            if (radius != -1.0) {
                visibility_region.radius = radius;
            }
            continue;
        }
        if (token == "[VERTICES]") {
            int n_vertices;
            ifs >> n_vertices;
            visibility_region.vertices.resize(n_vertices);
            for (int i = 0; i < n_vertices; ++i) {
                int vertex_id;
                ifs >> vertex_id;
                auto &v = visibility_region.vertices[vertex_id];
                ifs >> v.point.x;
                ifs >> v.point.y;
                ifs >> v.edge_flag;
                ifs >> v.vertex_flag;
            }
            continue;
        }
    }
    return true;
}

std::string data_loading::LoadVisibilityRegionSafely(
    const std::string &file,
    RadialVisibilityRegion &visibility_region
) {
    try {
        if (!LoadVisibilityRegion(file, visibility_region)) {
            return "File " + file + " could not be opened or found.";
        }
    } catch (const std::invalid_argument &e) {
        return "Error while loading " + file + ":\n" + e.what();
    }
    return "ok";
}
