/**
 * File:   visibility_region.h
 *
 * Date:   30.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_VISIBILITY_REGIONS_H_
#define TRIVIS_VISIBILITY_REGIONS_H_

#include <vector>

#include "trivis/geom/geom_types.h"

namespace trivis {

struct AbstractVisibilityRegionVertex {
    bool is_intersection = false;
    int id = -1;
    int id_b = -1;
    int id_c = -1;
    int id_d = -1;
};

bool operator==(const AbstractVisibilityRegionVertex &v1, const AbstractVisibilityRegionVertex &v2);

bool operator!=(const AbstractVisibilityRegionVertex &v1, const AbstractVisibilityRegionVertex &v2);

struct AbstractVisibilityRegionSegment {
    int id = -1;
    AbstractVisibilityRegionVertex v1;
    AbstractVisibilityRegionVertex v2;
};

bool operator==(const AbstractVisibilityRegionSegment &s1, const AbstractVisibilityRegionSegment &s2);

bool operator!=(const AbstractVisibilityRegionSegment &s1, const AbstractVisibilityRegionSegment &s2);

struct AbstractVisibilityRegion {
    int seed_id = -1;
    geom::FPoint seed;
    std::vector<AbstractVisibilityRegionSegment> segments;
};

bool operator==(const AbstractVisibilityRegion &r1, const AbstractVisibilityRegion &r2);

bool operator!=(const AbstractVisibilityRegion &r1, const AbstractVisibilityRegion &r2);

struct VisibilityRegionVertex {
    int vertex_flag;
    int edge_flag;
    geom::FPoint point;
};

struct RadialVisibilityRegion {
    double radius = -1;
    int seed_id = -1;
    geom::FPoint seed;
    std::vector<VisibilityRegionVertex> vertices;
};

bool IsValid(const RadialVisibilityRegion &reg);

void RemoveAntennas(
    RadialVisibilityRegion &res
);

double Area(
    const RadialVisibilityRegion &vis_reg
);

void RemoveShortEdges(
    double min_edge_length,
    RadialVisibilityRegion &res
);

RadialVisibilityRegion SampleArcEdges(
    const RadialVisibilityRegion &reg,
    double max_sample_beta
);

}

#endif //TRIVIS_VISIBILITY_REGIONS_H_
