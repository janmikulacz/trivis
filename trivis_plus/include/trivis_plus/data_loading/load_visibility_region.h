/**
 * File:   load_visibility_region.h
 *
 * Date:   27.04.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PLUS_DATA_LOADING_LOAD_VISIBILITY_REGION_H_
#define TRIVIS_PLUS_DATA_LOADING_LOAD_VISIBILITY_REGION_H_

#include <string>
#include <optional>

#include "trivis/visibility_regions.h"

namespace trivis_plus::data_loading {

bool LoadVisibilityRegion(
    const std::string &file,
    trivis::RadialVisibilityRegion &visibility_region
) noexcept(false);

std::string LoadVisibilityRegionSafely(
    const std::string &file,
    trivis::RadialVisibilityRegion &visibility_region
);

}

#endif //TRIVIS_PLUS_DATA_LOADING_LOAD_VISIBILITY_REGION_H_
