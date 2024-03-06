/**
 * File:   vis_regions.h
 *
 * Date:   30.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_VIS_REGIONS_H_
#define TRIVIS_VIS_REGIONS_H_

#include <vector>
#include <optional>

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
    std::optional<double> radius;
    std::optional<int> seed_id;
    geom::FPoint seed;
    std::vector<VisibilityRegionVertex> vertices;

    [[nodiscard]] bool IsValid() const;

    [[nodiscard]] double Area() const;

    void RemoveAntennas();

    void IntersectWithCircleCenteredAtSeed(std::optional<double> circle_radius = std::nullopt);

    void RemoveShortEdges(double min_edge_length);

    void SampleArcEdges(double max_sample_angle);

    [[nodiscard]] geom::FPolygon ToPolygon() const;
};

bool IsValid(const RadialVisibilityRegion &region);

double Area(const RadialVisibilityRegion &region);

void RemoveAntennas(RadialVisibilityRegion &region);

void IntersectWithCircleCenteredAtSeed(RadialVisibilityRegion &region, std::optional<double> radius = std::nullopt);

void RemoveShortEdges(RadialVisibilityRegion &region, double min_edge_length);

void SampleArcEdges(RadialVisibilityRegion &region, double max_sample_angle);

geom::FPolygon ToPolygon(const RadialVisibilityRegion &region);

inline bool operator==(const AbstractVisibilityRegionVertex &v1, const AbstractVisibilityRegionVertex &v2) {
    bool is_intersection = v1.is_intersection;
    if (v2.is_intersection != is_intersection) {
        return false; // The intersection flag must be the same.
    }
    if (!is_intersection) {
        return v1.id == v2.id; // If not intersection, then id_b, id_c, id_d does not matter.
    }
    return v1.id == v2.id && v1.id_b == v2.id_b && v1.id_c == v2.id_c && v1.id_d == v2.id_d;
}

inline bool operator!=(const AbstractVisibilityRegionVertex &v1, const AbstractVisibilityRegionVertex &v2) {
    return !(v1 == v2);
}

inline bool operator==(const AbstractVisibilityRegionSegment &s1, const AbstractVisibilityRegionSegment &s2) {
    return s1.id == s2.id && s1.v1 == s2.v1 && s1.v2 == s2.v2;
}

inline bool operator!=(const AbstractVisibilityRegionSegment &s1, const AbstractVisibilityRegionSegment &s2) {
    return !(s1 == s2);
}

inline bool operator==(const AbstractVisibilityRegion &r1, const AbstractVisibilityRegion &r2) {
    return r1.seed_id == r2.seed_id && r1.seed == r2.seed && r1.segments == r2.segments;
}

inline bool operator!=(const AbstractVisibilityRegion &r1, const AbstractVisibilityRegion &r2) {
    return !(r1 == r2);
}

inline bool RadialVisibilityRegion::IsValid() const {
    return trivis::IsValid(*this);
}

inline double RadialVisibilityRegion::Area() const {
    return trivis::Area(*this);
}

inline void RadialVisibilityRegion::RemoveAntennas() {
    trivis::RemoveAntennas(*this);
}

inline void RadialVisibilityRegion::IntersectWithCircleCenteredAtSeed(std::optional<double> circle_radius) {
    trivis::IntersectWithCircleCenteredAtSeed(*this, circle_radius);
}

inline void RadialVisibilityRegion::RemoveShortEdges(double min_edge_length) {
    trivis::RemoveShortEdges(*this, min_edge_length);
}

inline void RadialVisibilityRegion::SampleArcEdges(double max_sample_angle) {
    trivis::SampleArcEdges(*this, max_sample_angle);
}

inline geom::FPolygon RadialVisibilityRegion::ToPolygon() const {
    return trivis::ToPolygon(*this);
}

}

#endif //TRIVIS_VIS_REGIONS_H_
