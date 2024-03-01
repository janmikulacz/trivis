/**
 * File:   load_map.h
 *
 * Date:   29.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DATA_LOADING_LOAD_MAP_H_
#define TRIVIS_PLUS_DATA_LOADING_LOAD_MAP_H_

#include <string>
#include <optional>

#include "trivis/geom/poly_map.h"

namespace trivis_plus::data_loading {

/**
 *
 * Loads map from file in the following format:
 *
 *      [SCALE]
 *      <default_scale>
 *
 *      [BORDER]
 *      <x1> <y1>
 *      <x2> <y2>
 *      ...
 *      <xn> <yn>
 *
 *      [OBSTACLE]
 *      <x1> <y1>
 *      <x2> <y2>
 *      ...
 *      <xi> <yi>
 *
 *      ...
 *
 *      [OBSTACLE]
 *      <x1> <y1>
 *      <x2> <y2>
 *      ...
 *      <xj> <yj>
 *
 * @param file
 * @param border
 * @param holes
 * @param scale
 * @return
 */
bool LoadPolyMap(
    const std::string &file,
    trivis::geom::PolyMap &poly_map,
    std::optional<double> scale = std::nullopt
) noexcept(false);

std::string LoadPolyMapSafely(
    const std::string &file,
    trivis::geom::PolyMap &poly_map,
    std::optional<double> scale = std::nullopt
);

}

#endif //TRIVIS_PLUS_DATA_LOADING_LOAD_MAP_H_
