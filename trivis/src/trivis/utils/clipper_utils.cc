/**
 * File:   clipper_utils.cc
 *
 * Date:   28.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include "trivis/utils/clipper_utils.h"

using namespace trivis;
using namespace trivis::geom;
using namespace trivis::utils;

Clipper2Lib::PathD utils::Geom2Clipper(
    const FPolygon &polygon,
    double normalizer,
    const FPoint &subtrahend,
    double multiplier
) {
    Clipper2Lib::PathD ret;
    ret.reserve(polygon.size());
    for (const auto &p: polygon) {
        ret.push_back(Geom2Clipper(p, normalizer, subtrahend, multiplier));
    }
    return ret;
}

FPolygon utils::Clipper2Geom(
    const Clipper2Lib::PathD &polygon,
    double normalizer,
    const FPoint &subtrahend,
    double multiplier
) {
    FPolygon ret;
    ret.reserve(polygon.size());
    for (const auto &p: polygon) {
        ret.push_back(Clipper2Geom(p, normalizer, subtrahend, multiplier));
    }
    return ret;
}

Clipper2Lib::PathsD utils::Geom2Clipper(
    const FPolygons &polygons,
    double normalizer,
    const FPoint &subtrahend,
    double multiplier
) {
    Clipper2Lib::PathsD ret;
    ret.reserve(polygons.size());
    for (const auto &polygon: polygons) {
        ret.push_back(Geom2Clipper(polygon, normalizer, subtrahend, multiplier));
    }
    return ret;
}

FPolygons utils::Clipper2Geom(
    const Clipper2Lib::PathsD &polygons,
    double normalizer,
    const FPoint &subtrahend,
    double multiplier
) {
    FPolygons ret;
    ret.reserve(polygons.size());
    for (const auto &polygon: polygons) {
        ret.push_back(Clipper2Geom(polygon, normalizer, subtrahend, multiplier));
    }
    return ret;
}