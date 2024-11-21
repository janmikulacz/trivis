/**
 * File:   utils.cc
 *
 * Date:   08.11.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis_pathfinder/utils/utils.h"

using namespace trivis_pathfinder;
using namespace trivis_pathfinder::utils;

void utils::ConvertPathToPoints(
    const std::vector<int> &cities_id_path,
    const trivis::Trivis &vis,
    trivis_pathfinder::TrivisPathfinder &path_finder,
    std::vector<int> &cities_reflex_id_path,
    std::vector<bool> &cities_reflex_id_path_is_city,
    trivis::geom::FPoints &cities_reflex_point_path
) {
    cities_reflex_id_path.clear();
    cities_reflex_id_path_is_city.clear();
    cities_reflex_point_path.clear();
    int n_path = static_cast<int>(cities_id_path.size());
    for (int i_prev = 0, i = 1; i < n_path; i_prev = i++) {
        int city_id_i_prev = cities_id_path[i_prev];
        int city_id_i = cities_id_path[i];
        trivis::geom::FPoints point_path;
        std::vector<int> id_path;
        path_finder.ShortestPathCities(vis, city_id_i_prev, city_id_i, &point_path, &id_path);
        if (i + 1 == n_path) {
            id_path.push_back(path_finder.n_reflex() + city_id_i);
        } else {
            point_path.pop_back();
        }
        cities_reflex_id_path.push_back(path_finder.n_reflex() + city_id_i_prev);
        cities_reflex_id_path = trivis::utils::Concatenate(cities_reflex_id_path, id_path);
        cities_reflex_id_path_is_city.push_back(true);
        cities_reflex_id_path_is_city = trivis::utils::Concatenate(cities_reflex_id_path_is_city, std::vector<bool>(id_path.size(), false));
        if (i + 1 == n_path) {
            cities_reflex_id_path_is_city.back() = true;
        }
        cities_reflex_point_path = trivis::utils::Concatenate(cities_reflex_point_path, point_path);
    }
}