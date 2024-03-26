/**
 * File:   load_map.cc
 *
 * Date:   29.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_plus/data_loading/load_map.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace trivis_plus;
using namespace trivis_plus::data_loading;

using namespace trivis;
using namespace trivis::geom;

bool data_loading::SavePolyMap(
    const PolyMap &poly_map,
    const std::string &file,
    std::stringstream *info
) {
    std::ofstream ofs(file);
    if (!ofs.is_open()) {
        if (info) *info << "Cannot open file " << file << " for writing.\n";
        return false;
    }
    ofs << "[SCALE]\n";
    ofs << "1.0\n";
    ofs << "\n[BORDER]\n";
    for (const auto &p: poly_map.border()) {
        ofs << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << p.x;
        ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << p.y;
        ofs << "\n";
    }
    for (const auto &hole: poly_map.holes()) {
        ofs << "\n[OBSTACLE]\n";
        for (const auto &p: hole) {
            ofs << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << p.x;
            ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << p.y;
            ofs << "\n";
        }
    }
    return true;
}

void TryLoadPolyMapNewFormat(
    std::ifstream &ifs,
    PolyMap &poly_map
) noexcept(false) {
    std::string token;
    FPoints points;
    FPolygon border;
    FPolygons holes;
    while (true) {
        ifs >> token;
        if (ifs.eof()) {
            break;
        }
        if (token == "[NAME]") {
            // can be ignored
            continue;
        }
        if (token == "[LIMITS]") {
            // can be ignored
            continue;
        }
        if (token == "[POINTS]") {
            int n_points;
            ifs >> n_points;
            points.reserve(n_points);
            int id;
            double x, y;
            for (int i = 0; i < n_points; ++i) {
                ifs >> id;
                ifs >> x;
                ifs >> y;
                points.emplace_back(x, y);
            }
            continue;
        }
        if (token == "[BORDER]") {
            int n_border;
            ifs >> n_border;
            border.reserve(n_border);
            int id;
            int point_id;
            for (int i = 0; i < n_border; ++i) {
                ifs >> id;
                ifs >> point_id;
                border.push_back(points[point_id]);
            }
            continue;
        }
        if (token == "[NUM_HOLES]") {
            int n_holes;
            ifs >> n_holes;
            holes.resize(n_holes);
            continue;
        }
        if (token == "[HOLE]") {
            int hole_id, n_hole;
            ifs >> hole_id;
            ifs >> n_hole;
            auto &hole = holes[hole_id];
            hole.reserve(n_hole);
            int id;
            int point_id;
            for (int i = 0; i < n_hole; ++i) {
                ifs >> id;
                ifs >> point_id;
                hole.push_back(points[point_id]);
            }
            continue;
        }
    }
    // Save to output.
    poly_map = PolyMap(border, holes);
}

void TryLoadPolyMap(
    std::ifstream &ifs,
    PolyMap &poly_map,
    std::optional<double> scale
) noexcept(false) {
    std::string token;
    FPolygon curr_polygon;
    bool processing_border = false;
    bool processing_obstacle = false;
    bool processing_map_points = false;
    bool format_map_points = false;
    Points<double> points;
    FPolygon border;
    FPolygons holes;
    double d_max = std::numeric_limits<double>::max();
    double x_min = -d_max, x_max = d_max, y_min = -d_max, y_max = d_max;
    while (true) {
        ifs >> token;
        if (token == "[NAME]") {
            TryLoadPolyMapNewFormat(ifs, poly_map);
            return;
        }
        if (token == "[INFO]") {
            continue;
        }
        if (token == "FORMAT=MAP_POINTS") {
            format_map_points = true;
            if (!scale || scale < 0.0) {
                scale = 1.0;
            }
            continue;
        }
        if (!format_map_points) { // FORMAT: SIMPLE
            // Save the border flag.
            if (token == "[BORDER]") {
                processing_border = true;
            }
            // If the next border/obstacle polygon or the end of the file is detected, ...
            // ... then save the current polygon to the map if it is non-empty.
            if (token == "[BORDER]" || token == "[OBSTACLE]" || ifs.eof()) {
                if (!curr_polygon.empty()) {
                    // Check polygon orientations.
                    if (processing_border) {
                        // Border should return true.
                        if (!OrientationCounterClockwise(curr_polygon)) {
                            // Fix the orientation if wrong.
                            ChangeOrientation(curr_polygon);
                        }
                        processing_border = false;
                        border = curr_polygon;
                        ComputeLimits(border, x_min, x_max, y_min, y_max);
                    } else {
                        // Obstacle should return false.
                        if (OrientationCounterClockwise(curr_polygon)) {
                            // Fix the orientation if wrong.
                            ChangeOrientation(curr_polygon);
                        }
                        holes.push_back(curr_polygon);
                    }
                    // Add the polygon to the map.
                    curr_polygon.clear();
                    if (ifs.eof()) break;
                }
            } else if (token == "[SCALE]") {
                ifs >> token;
                if (!scale || scale < 0.0) { // Load scale only if it is not given.
                    scale = std::stod(token);
                }
            } else if (!ifs.eof()) { // Assuming line with two double coordinates.
                double x, y;
                x = std::stod(token) * *scale;
                ifs >> token;
                y = std::stod(token) * *scale;
                curr_polygon.emplace_back(std::min(std::max(x_min, x), x_max), std::min(std::max(y_min, y), y_max));
            }
        } else { // FORMAT: MAP POINTS

            if (token == "[MAP_BORDER]" || token == "[MAP_OBSTACLE]" || ifs.eof()) {
                if (!curr_polygon.empty()) {
                    // Check polygon orientations.
                    if (processing_border) {
                        // Border should return true.
                        if (!OrientationCounterClockwise(curr_polygon)) {
                            // Fix the orientation if wrong.
                            ChangeOrientation(curr_polygon);
                        }
                        processing_border = false;
                        border = curr_polygon;
                        ComputeLimits(border, x_min, x_max, y_min, y_max);
                    } else {
                        // Obstacle should return false.
                        if (OrientationCounterClockwise(curr_polygon)) {
                            // Fix the orientation if wrong.
                            ChangeOrientation(curr_polygon);
                        }
                        holes.push_back(curr_polygon);
                    }
                    curr_polygon.clear();
                }
                if (ifs.eof()) {
                    break;
                }
            }

            if (token == "[MAP_POINTS]") {
                processing_map_points = true;
                processing_border = false;
                processing_obstacle = false;
            } else if (token == "[MAP_BORDER]") {
                processing_map_points = false;
                processing_border = true;
                processing_obstacle = false;
            } else if (token == "[MAP_OBSTACLE]") {
                processing_map_points = false;
                processing_border = false;
                processing_obstacle = true;
            } else if (token == "[MAP_CONVEX_REGION]" ||
                token == "[MESH_NODES]" ||
                token == "[MESH_TRIANGLES]" ||
                token == "[MESH_BOUNDARY_NODES]" ||
                token == "[TRIMESH_CONVEX_REGION]" ||
                token == "[TRIMESH_CONVEX_REGION_TRIANGLES]") {
                processing_map_points = false;
                processing_border = false;
                processing_obstacle = false;
            } else if (processing_map_points) {
                double x, y;
                ifs >> x;
                ifs >> y;
                points.emplace_back(std::min(std::max(x_min, x * *scale), x_max), std::min(std::max(y_min, y * *scale), y_max));
            } else if (processing_border || processing_obstacle) {
                int idx = std::stoi(token);
                curr_polygon.push_back(points[idx]);
            }
        }
    }
    poly_map = PolyMap(border, holes);
}

std::optional<trivis::geom::PolyMap> data_loading::LoadPolyMap(
    const std::string &file,
    std::optional<double> rescale,
    std::stringstream *info
) {
    std::ifstream ifs(file.c_str());
    if (ifs.fail()) {
        if (info) *info << "Cannot open file " << file << " for reading.\n";
        return std::nullopt;
    }
    trivis::geom::PolyMap poly_map;
    try {
        TryLoadPolyMap(ifs, poly_map, rescale);
    } catch (const std::exception &e) {
        if (info) *info << "Error while loading " << file << ":\n" << e.what();
        return std::nullopt;
    }
    return poly_map;
}
