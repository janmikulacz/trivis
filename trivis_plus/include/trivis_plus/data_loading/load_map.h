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
#include <sstream>
#include <optional>

#include "trivis/geom/poly_map.h"

namespace trivis_plus::data_loading {

bool SavePolyMap(
    const trivis::geom::PolyMap &poly_map,
    const std::string &file,
    std::stringstream *info = nullptr
);

std::optional<trivis::geom::PolyMap> LoadPolyMap(
    const std::string &file,
    std::optional<double> rescale = std::nullopt,
    std::stringstream *info = nullptr
);

}

#endif //TRIVIS_PLUS_DATA_LOADING_LOAD_MAP_H_
