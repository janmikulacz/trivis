/**
 * File:   random_points.cc
 *
 * Date:   06.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/utils/random_points.h"

#include "trivis/geom/generic_geom_utils.h"

#include <cassert>

using namespace trivis;
using namespace trivis::geom;
using namespace trivis::utils;

FPoint utils::UniformRandomPointInTriangle(
    const FPolygon &triangle,
    std::mt19937 &rng
) {
    auto distribution_01 = std::uniform_real_distribution<double>(0.0, 1.0);
    auto r1_sqrt = std::sqrt(distribution_01(rng));
    auto r2 = distribution_01(rng);
    assert(triangle.size() >= 3);
    auto x = (1 - r1_sqrt) * triangle[0].x + r1_sqrt * (1 - r2) * triangle[1].x + r2 * r1_sqrt * triangle[2].x;
    auto y = (1 - r1_sqrt) * triangle[0].y + r1_sqrt * (1 - r2) * triangle[1].y + r2 * r1_sqrt * triangle[2].y;
    return geom::MakePoint(x, y);
}

FPoint utils::UniformRandomPointInRandomTriangle(
    const FPolygons &triangles,
    const std::vector<double> &accum_areas,
    std::mt19937 &rng,
    int *triangle_id
) {
    auto r = std::uniform_real_distribution<double>(0.0, accum_areas.back())(rng);
    int r_idx = 0;
    for (; r_idx < triangles.size(); ++r_idx) {
        if (accum_areas[r_idx] >= r) {
            break;
        }
    }
    const auto &triangle = triangles[r_idx];
    auto distribution_01 = std::uniform_real_distribution<double>(0.0, 1.0);
    auto r1_sqrt = std::sqrt(distribution_01(rng));
    auto r2 = distribution_01(rng);
    assert(triangle.size() >= 3);
    auto x = (1 - r1_sqrt) * triangle[0].x + r1_sqrt * (1 - r2) * triangle[1].x + r2 * r1_sqrt * triangle[2].x;
    auto y = (1 - r1_sqrt) * triangle[0].y + r1_sqrt * (1 - r2) * triangle[1].y + r2 * r1_sqrt * triangle[2].y;
    if (triangle_id) *triangle_id = r_idx;
    return geom::MakePoint(x, y);
}

FPoint utils::UniformRandomPointInRandomTriangle(
    const FPolygons &triangles,
    std::vector<double> &accum_areas,
    std::mt19937 &rng,
    int *triangle_id
) {
    if (!accum_areas.empty()) {
        assert(triangles.size() == accum_areas.size());
    } else {
        for (const auto &triangle: triangles) {
            accum_areas.push_back((accum_areas.empty() ? 0.0 : accum_areas.back()) + Area(triangle));
        }
    }
    return UniformRandomPointInRandomTriangle(triangles, const_cast<const std::vector<double> &>(accum_areas), rng, triangle_id);
}
