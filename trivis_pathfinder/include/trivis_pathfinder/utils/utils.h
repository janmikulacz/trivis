/**
 * File:   utils.h
 *
 * Date:   08.11.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PATHFINDER_UTILS_H_
#define TRIVIS_PATHFINDER_UTILS_H_

#include "trivis_pathfinder/trivis_pathfinder.h"

namespace trivis_pathfinder::utils {

void ConvertPathToPoints(
    const std::vector<int> &cities_id_path,
    const trivis::Trivis &vis,
    TrivisPathfinder &path_finder,
    std::vector<int> &cities_reflex_id_path,
    std::vector<bool> &cities_reflex_id_path_is_city,
    trivis::geom::FPoints &cities_reflex_point_path
);

}

#endif //TRIVIS_PATHFINDER_UTILS_H_
