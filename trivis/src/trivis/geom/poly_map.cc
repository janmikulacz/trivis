/**
 * File:   poly_map.cc
 *
 * Date:   23.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/geom/poly_map.h"

#include <map>

#include "trivis/geom/robust_geometry.h"

using namespace trivis;
using namespace trivis::geom;

void PolyMap::RemoveDuplicatePoints() {
    _border = geom::RemoveDuplicatePoints(_border);
    for (auto &hole: _holes) {
        hole = geom::RemoveDuplicatePoints(hole);
    }
}

void PolyMap::RemoveCollinearPoints() {
    _border = geom::RemoveCollinearPoints(_border);
    for (auto &hole: _holes) {
        hole = geom::RemoveCollinearPoints(hole);
    }
}

void PolyMap::ShiftToOrigin() {
    for (auto &p: _border) {
        p.x -= _limits.x_min;
        p.y -= _limits.y_min;
    }
    for (auto &polygon: _holes) {
        for (auto &p: polygon) {
            p.x -= _limits.x_min;
            p.y -= _limits.y_min;
        }
    }
    _limits.x_max -= _limits.x_min;
    _limits.y_max -= _limits.y_min;
    _limits.x_min = 0.0;
    _limits.y_min = 0.0;
}
FPoints PolyMap::ToPoints() const {
    FPoints ret;
    for (const auto &p: _border) {
        ret.push_back(p);
    }
    for (const auto &hole: _holes) {
        for (const auto &p: hole) {
            ret.push_back(p);
        }
    }
    return ret;
}

FPolygons SimplifyWeaklySimplePolygon(
    const FPolygon &polygon
) {
    std::map<std::pair<double, double>, int> unique_points;
    for (int i = 0; i < polygon.size(); ++i) {
        const auto &p = polygon[i];
        const auto p_pair = std::make_pair(p.x, p.y);
        auto it = unique_points.find(p_pair);
        if (it != unique_points.end()) {
            FPolygon remaining_polygon(polygon.begin(), polygon.begin() + it->second);
            FPolygon simple_polygon(polygon.begin() + it->second, polygon.begin() + i);
            remaining_polygon.insert(remaining_polygon.end(), polygon.begin() + i, polygon.end());
            FPolygons ret = SimplifyWeaklySimplePolygon(remaining_polygon);
            ret.push_back(simple_polygon);
            return ret;
        }
        unique_points.emplace(p_pair, i);
    }
    return {polygon};
}

void PolyMap::SimplifyWeaklySimplePolygons() {
    auto ret = SimplifyWeaklySimplePolygon(_border);
    _border = std::move(ret[0]);
    FPolygons new_holes;
    for (int k = 1; k < ret.size(); ++k) {
        if (OrientationClockwise(ret[k])) {
            new_holes.push_back(std::move(ret[k]));
        }
    }
    for (const auto &hole: _holes) {
        ret = SimplifyWeaklySimplePolygon(hole);
        for (auto &new_hole: ret) {
            if (OrientationClockwise(new_hole)) {
                new_holes.push_back(std::move(new_hole));
            }
        }
    }
    _holes = std::move(new_holes);
}
