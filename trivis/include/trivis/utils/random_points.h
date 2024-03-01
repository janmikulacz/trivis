/**
 * File:   random_points.h
 *
 * Date:   06.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_UTILS_RANDOM_POINTS_H_
#define TRIVIS_UTILS_RANDOM_POINTS_H_

#include <random>

#include "trivis/geom/geom_types.h"

namespace trivis::utils {

[[nodiscard]] geom::FPoint UniformRandomPointInTriangle(
    const geom::FPolygon &triangle,
    std::mt19937 &rng
);

[[nodiscard]] geom::FPoint UniformRandomPointInRandomTriangle(
    const geom::FPolygons &triangles,
    const std::vector<double> &accum_areas,
    std::mt19937 &rng,
    int *triangle_id = nullptr
);

[[nodiscard]] geom::FPoint UniformRandomPointInRandomTriangle(
    const geom::FPolygons &triangles,
    std::vector<double> &accum_areas,
    std::mt19937 &rng,
    int *triangle_id = nullptr
);

[[nodiscard]] inline geom::FPoint UniformRandomPointInRandomTriangle(
    const geom::FPolygons &triangles,
    std::mt19937 &rng,
    int *triangle_id = nullptr
) {
    std::vector<double> accum_areas;
    return UniformRandomPointInRandomTriangle(triangles, accum_areas, rng, triangle_id);
}

}

#endif //TRIVIS_UTILS_RANDOM_POINTS_H_
