/**
 * File:   clipper_geom.h
 *
 * Date:   31.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_UTILS_CLIPPER_GEOM_H_
#define TRIVIS_UTILS_CLIPPER_GEOM_H_

#include "trivis/geom/poly_map.h"
#include "trivis/vis_regions.h"
#include "trivis/utils/clipper_utils.h"

namespace trivis::utils {

Clipper2Lib::Paths64 ToClipper(
    const geom::PolyMap &map
);

int64_t ToClipper(
    double x,
    const geom::FLimits &lim
);

Clipper2Lib::Paths64 ToClipper(
    const geom::FPolygons &polygons,
    const geom::FLimits &lim
);

Clipper2Lib::Path64 ToClipper(
    const geom::FPolygon &polygon,
    const geom::FLimits &lim
);

Clipper2Lib::Path64 ToClipper(
    const RadialVisibilityRegion &vis_reg,
    const geom::FLimits &lim
);

double FromClipper(
    int64_t x,
    const geom::FLimits &lim
);

geom::FPolygon FromClipper(
    const Clipper2Lib::Path64 &polygon,
    const geom::FLimits &lim
);

geom::FPolygons FromClipper(
    const Clipper2Lib::Paths64 &polygons,
    const geom::FLimits &lim
);

void ClipOff(
    const Clipper2Lib::Path64 &region,
    Clipper2Lib::Paths64 &map
);

void ClipOff(
    const Clipper2Lib::Paths64 &regions,
    Clipper2Lib::Paths64 &map
);

// void ClipOff(
//     const std::vector<VisibilityRegionApproximations> &coverage,
//     Clipper2Lib::Paths64 &map
// );

void ClipOff(
    const Clipper2Lib::Path64 &region,
    std::vector<Clipper2Lib::Paths64> &map
);

void ClipOff(
    const Clipper2Lib::Paths64 &regions,
    std::vector<Clipper2Lib::Paths64> &map
);

// void ClipOff(
//     const std::vector<VisibilityRegionApproximations> &coverage,
//     std::vector<Clipper2Lib::Paths64> &map
// );

Clipper2Lib::Paths64 Intersection(
    const Clipper2Lib::Path64 &region,
    const Clipper2Lib::Paths64 &map
);

double ClipperArea(
    const Clipper2Lib::Paths64 &regions
);

double ClipperArea(
    const std::vector<Clipper2Lib::Paths64> &regions
);

// Clipper2Lib::Paths64 Union(
//     const std::vector<VisibilityRegionApproximations> &subjects,
//     int cluster_size = 100
// );

// Clipper2Lib::Paths64 Difference(
//     const Clipper2Lib::Paths64 &subject,
//     const VisibilityRegionApproximations &clip
// );

// Clipper2Lib::Paths64 Difference(
//     const Clipper2Lib::Paths64 &subject,
//     const std::vector<VisibilityRegionApproximations> &clips,
//     int cluster_size = 100
// );

}

#endif //TRIVIS_UTILS_CLIPPER_GEOM_H_
