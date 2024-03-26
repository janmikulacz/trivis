/**
 * File:   load_vis_region.h
 *
 * Date:   27.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DATA_LOADING_LOAD_VIS_REGION_H_
#define TRIVIS_PLUS_DATA_LOADING_LOAD_VIS_REGION_H_

#include <string>
#include <sstream>
#include <optional>

#include "trivis/vis_regions.h"

namespace trivis_plus::data_loading {

bool SaveVisRegion(
    const trivis::RadialVisibilityRegion &region,
    const std::string &file,
    std::stringstream *info = nullptr
);

std::optional<trivis::RadialVisibilityRegion> LoadVisRegion(
    const std::string &file,
    std::stringstream *info = nullptr
);

}

#endif //TRIVIS_PLUS_DATA_LOADING_LOAD_VIS_REGION_H_
