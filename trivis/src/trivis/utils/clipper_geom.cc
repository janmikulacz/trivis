/**
 * File:   clipper_geom.cc
 *
 * Date:   31.03.2023
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/utils/clipper_geom.h"

using namespace trivis;
using namespace trivis::geom;
using namespace trivis::utils;

Clipper2Lib::Paths64 utils::ToClipper(
    const PolyMap &map
) {
    Clipper2Lib::Paths64 ret;
    const auto &lim = map.limits();
    double normalizer = std::max(lim.x_max - lim.x_min, lim.y_max - lim.y_min);
    ret.reserve(1 + map.holes().size());
    ret.push_back(Clipper2Lib::PathDToPath64(utils::Geom2Clipper(map.border(), normalizer, {lim.x_min, lim.y_min})));
    for (const auto &hole: map.holes()) {
        ret.push_back(Clipper2Lib::PathDToPath64(utils::Geom2Clipper(hole, normalizer, {lim.x_min, lim.y_min})));
    }
    return ret;
}

int64_t utils::ToClipper(double x, const FLimits &lim) {
    double normalizer = std::max(lim.x_max - lim.x_min, lim.y_max - lim.y_min);
    return static_cast<int64_t>(utils::Number2Clipper(x, normalizer));
}

Clipper2Lib::Paths64 utils::ToClipper(const FPolygons &polygons, const FLimits &lim) {
    Clipper2Lib::Paths64 ret;
    ret.reserve(polygons.size());
    for (const auto &polygon: polygons) {
        ret.push_back(ToClipper(polygon, lim));
    }
    return ret;
}

Clipper2Lib::Path64 utils::ToClipper(const FPolygon &polygon, const FLimits &lim) {
    double normalizer = std::max(lim.x_max - lim.x_min, lim.y_max - lim.y_min);
    return Clipper2Lib::PathDToPath64(utils::Geom2Clipper(polygon, normalizer, {lim.x_min, lim.y_min}));
}

Clipper2Lib::Path64 utils::ToClipper(
    const RadialVisibilityRegion &vis_reg,
    const FLimits &lim
) {
    FPolygon polygon;
    polygon.reserve(vis_reg.vertices.size());
    for (const auto &v: vis_reg.vertices) {
        polygon.push_back(v.point);
    }
    double normalizer = std::max(lim.x_max - lim.x_min, lim.y_max - lim.y_min);
    return Clipper2Lib::PathDToPath64(utils::Geom2Clipper(polygon, normalizer, {lim.x_min, lim.y_min}));
}

double utils::FromClipper(
    int64_t x,
    const geom::FLimits &lim
) {
    double normalizer = std::max(lim.x_max - lim.x_min, lim.y_max - lim.y_min);
    return utils::NumberFromClipper(static_cast<double>(x), normalizer);
}

FPolygon utils::FromClipper(
    const Clipper2Lib::Path64 &polygon,
    const FLimits &lim
) {
    double normalizer = std::max(lim.x_max - lim.x_min, lim.y_max - lim.y_min);
    return utils::Clipper2Geom(Clipper2Lib::Path64ToPathD(polygon), normalizer, {lim.x_min, lim.y_min});
}

FPolygons utils::FromClipper(
    const Clipper2Lib::Paths64 &polygons,
    const FLimits &lim
) {
    FPolygons ret;
    ret.reserve(polygons.size());
    for (const auto &polygon: polygons) {
        ret.push_back(FromClipper(polygon, lim));
    }
    return ret;
}

void utils::ClipOff(
    const Clipper2Lib::Path64 &region,
    Clipper2Lib::Paths64 &map
) {
    Clipper2Lib::Clipper64 clipper;
    clipper.AddSubject(map);
    clipper.AddClip({region});
    clipper.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, map);
    clipper.Clear();
}

void utils::ClipOff(
    const Clipper2Lib::Paths64 &regions,
    Clipper2Lib::Paths64 &map
) {
    Clipper2Lib::Clipper64 clipper;
    clipper.AddSubject(map);
    clipper.AddClip(regions);
    clipper.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, map);
    clipper.Clear();
}

void utils::ClipOff(
    const Clipper2Lib::Path64 &region,
    std::vector<Clipper2Lib::Paths64> &map
) {
    Clipper2Lib::Clipper64 clipper;
    for (const auto &reg: map) {
        clipper.AddSubject(reg);
    }
    clipper.AddClip({region});
    Clipper2Lib::PolyTree64 sol;
    clipper.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, sol);
    clipper.Clear();
    map.clear();
    for (int i = 0; i < sol.Count(); ++i) {
        const auto *curr = sol[i];
        const auto *parent = curr->Parent();
        Clipper2Lib::Paths64 aux = {curr->Polygon()};
        for (auto it = curr->begin(); it != curr->end(); ++it) {
            if (curr->Parent() == parent) {
                break;
            }
            aux.push_back(curr->Polygon());
        }
        std::vector<Clipper2Lib::Paths64> polygons;
        Clipper2Lib::Paths64 holes;
        for (const auto &reg: aux) {
            if (Clipper2Lib::IsPositive(reg)) {
                polygons.push_back({reg});
            } else {
                holes.push_back(reg);
            }
        }
        for (const auto &hole: holes) {
            for (auto &reg: polygons) {
                if (Clipper2Lib::PointInPolygon(hole.front(), reg.front()) == Clipper2Lib::PointInPolygonResult::IsInside) {
                    reg.push_back(hole);
                    break;
                }
            }
        }
        for (const auto &polygon: polygons) {
            map.push_back(polygon);
        }
    }
}

void utils::ClipOff(
    const Clipper2Lib::Paths64 &regions,
    std::vector<Clipper2Lib::Paths64> &map
) {
    Clipper2Lib::Clipper64 clipper;
    for (const auto &reg: map) {
        clipper.AddSubject(reg);
    }
    clipper.AddClip(regions);
    Clipper2Lib::PolyTree64 sol;
    clipper.Execute(Clipper2Lib::ClipType::Difference, Clipper2Lib::FillRule::NonZero, sol);
    clipper.Clear();
    map.clear();
    for (int i = 0; i < sol.Count(); ++i) {
        const auto *curr = sol[i];
        const auto *parent = curr->Parent();
        Clipper2Lib::Paths64 aux = {curr->Polygon()};
        for (auto it = curr->begin(); it != curr->end(); ++it) {
            if (curr->Parent() == parent) {
                break;
            }
            aux.push_back(curr->Polygon());
        }
        std::vector<Clipper2Lib::Paths64> polygons;
        Clipper2Lib::Paths64 holes;
        for (const auto &reg: aux) {
            if (Clipper2Lib::IsPositive(reg)) {
                polygons.push_back({reg});
            } else {
                holes.push_back(reg);
            }
        }
        for (const auto &hole: holes) {
            for (auto &reg: polygons) {
                if (Clipper2Lib::PointInPolygon(hole.front(), reg.front()) == Clipper2Lib::PointInPolygonResult::IsInside) {
                    reg.push_back(hole);
                    break;
                }
            }
        }
        for (const auto &polygon: polygons) {
            map.push_back(polygon);
        }
    }
}

Clipper2Lib::Paths64 utils::Intersection(
    const Clipper2Lib::Path64 &region,
    const Clipper2Lib::Paths64 &map
) {
    Clipper2Lib::Paths64 ret;
    Clipper2Lib::Clipper64 clipper;
    clipper.AddSubject({region});
    clipper.AddClip(map);
    clipper.Execute(Clipper2Lib::ClipType::Intersection, Clipper2Lib::FillRule::NonZero, ret);
    clipper.Clear();
    return ret;
}

double utils::ClipperArea(const Clipper2Lib::Paths64 &regions) {
    double area = 0.0;
    for (const auto &reg: regions) {
        area += Clipper2Lib::Area(reg);
    }
    return area;
}

double utils::ClipperArea(const std::vector<Clipper2Lib::Paths64> &regions) {
    double area = 0.0;
    for (const auto &reg: regions) {
        area += ClipperArea(reg);
    }
    return area;
}
