/**
 * File:   load_vis_region.cc
 *
 * Date:   27.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/data_loading/load_vis_region.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace trivis_plus;
using namespace trivis_plus::data_loading;

using namespace trivis;
using namespace trivis::geom;

bool data_loading::SaveVisRegion(
    const trivis::RadialVisibilityRegion &region,
    const std::string &file,
    std::stringstream *info
) {
    std::ofstream ofs(file);
    if (!ofs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for writing.\n";
        return false;
    }
    ofs << "[SEED]\n";
    ofs << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << region.seed.x;
    ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << region.seed.y;
    ofs << " " << (region.seed_id ? region.seed_id.value() : -1);
    ofs << "\n";
    ofs << "\n[RADIUS]\n";
    ofs << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << (region.radius ? region.radius.value() : -1.0);
    ofs << "\n";
    ofs << "\n[VERTICES]\n";
    ofs << region.vertices.size() << "\n";
    for (int i = 0; i < region.vertices.size(); ++i) {
        const auto &v = region.vertices[i];
        ofs << i << " ";
        ofs << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << v.point.x;
        ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << v.point.y;
        ofs << " " << v.edge_flag;
        ofs << " " << v.vertex_flag;
        ofs << "\n";
    }
    return true;
}

void TryLoadVisRegion(
    std::ifstream &ifs,
    RadialVisibilityRegion &region
) noexcept(false) {
    std::string token;
    while (true) {
        ifs >> token;
        if (ifs.eof()) {
            break;
        }
        if (token == "[SEED]") {
            ifs >> region.seed.x;
            ifs >> region.seed.y;
            int seed_id;
            ifs >> seed_id;
            if (seed_id != -1) {
                region.seed_id = seed_id;
            }
            continue;
        }
        if (token == "[RADIUS]") {
            double radius;
            ifs >> radius;
            if (radius != -1.0) {
                region.radius = radius;
            }
            continue;
        }
        if (token == "[VERTICES]") {
            int n_vertices;
            ifs >> n_vertices;
            region.vertices.resize(n_vertices);
            for (int i = 0; i < n_vertices; ++i) {
                int vertex_id;
                ifs >> vertex_id;
                auto &v = region.vertices[vertex_id];
                ifs >> v.point.x;
                ifs >> v.point.y;
                ifs >> v.edge_flag;
                ifs >> v.vertex_flag;
            }
            continue;
        }
    }
}

std::optional<trivis::RadialVisibilityRegion> data_loading::LoadVisRegion(
    const std::string &file,
    std::stringstream *info
) {
    std::ifstream ifs(file);
    if (!ifs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for reading.\n";
        return std::nullopt;
    }
    RadialVisibilityRegion region;
    try {
        TryLoadVisRegion(ifs, region);
    } catch (const std::exception &e) {
        if (info) *info << "Error while loading " << file << ":\n" << e.what();
        return std::nullopt;
    }
    return region;
}
