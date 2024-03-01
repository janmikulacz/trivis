/**
 * File:    utils.h
 *
 * Date:    26.08.2020
 * Author:  Jan Mikula
 * E-mail:  jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_GENERIC_GEOM_UTILS_H_
#define TRIVIS_GENERIC_GEOM_UTILS_H_

#include "trivis/geom/generic_geom_types.h"

#include <limits>
#include <algorithm>
#include <type_traits>

#include <iostream>

namespace trivis::geom {

template<typename T>
[[nodiscard]] Points<T> PositionsToPoints(
    const Positions<T> &positions
) {
    Points<T> points;
    points.reserve(positions.size());
    for (const auto &pos: positions) points.push_back(pos.ToPoint());
    return points;
}

template<typename T>
[[nodiscard]] Points<T> PolygonsToPoints(
    const Polygons<T> &polygons
) {
    Points<T> points;
    for (const auto &polygon: polygons) points.insert(points.end(), polygon.begin(), polygon.end());
    return points;
}

template<typename T>
[[nodiscard]] Polygon<T> MakeRegularPolygon(
    const Point<T> &center,
    T radius,
    unsigned n
) {
    Polygon<T> poly;
    double angle = 2.0 * M_PI / n;
    double alpha = angle / 2.0;
    for (int i = 0; i < n; ++i) {
        alpha += angle;
        poly.emplace_back(static_cast<T>(center.x + radius * sin(alpha)),
                          static_cast<T>(center.y + radius * cos(alpha)));
    }
    return poly;
}

template<typename T>
void ComputeLimits(
    const Polygon<T> &polygon,
    T &x_min,
    T &x_max,
    T &y_min,
    T &y_max
) {
    x_min = std::numeric_limits<T>::max();
    y_min = std::numeric_limits<T>::max();
    x_max = std::numeric_limits<T>::lowest();
    y_max = std::numeric_limits<T>::lowest();
    for (const auto &p: polygon) {
        if (p.x < x_min) x_min = p.x;
        if (p.y < y_min) y_min = p.y;
        if (x_max < p.x) x_max = p.x;
        if (y_max < p.y) y_max = p.y;
    }
}

template<typename T>
void ComputeLimits(
    const Polygons<T> &polygons,
    T &x_min,
    T &x_max,
    T &y_min,
    T &y_max
) {
    x_min = std::numeric_limits<T>::max();
    y_min = std::numeric_limits<T>::max();
    x_max = std::numeric_limits<T>::lowest();
    y_max = std::numeric_limits<T>::lowest();
    for (const auto &poly: polygons) {
        for (const auto &p: poly) {
            if (p.x < x_min) x_min = p.x;
            if (p.y < y_min) y_min = p.y;
            if (x_max < p.x) x_max = p.x;
            if (y_max < p.y) y_max = p.y;
        }
    }
}

template<typename T, typename TT = double>
[[nodiscard]] TT Area2(
    const Polygon<T> &polygon
) {
    int n = static_cast<int>(polygon.size());
    if (n < 3) return 0.0;
    TT area = 0;
    for (int i = 0, j = n - 1; i < n; ++i) {
        area += static_cast<double>(polygon[i].x + polygon[j].x) * static_cast<double>(polygon[i].y - polygon[j].y);
        j = i;
    }
    return area;
}

template<typename T, typename TT = double>
[[nodiscard]] inline TT Area(
    const Polygon<T> &polygon
) {
    return Area2<TT>(polygon) / TT(2);
}

template<typename T, typename TT = double>
[[nodiscard]] inline TT Area(
    const Polygons<T> &polygons
) {
    TT ret = TT(0);
    for (int i = 0; i < polygons.size(); ++i) {
        ret += Area<TT>(polygons[i]);
    }
    return ret;
}

/**
 *
 * Based on:
 * https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
 *
 * @tparam T
 * @param polygon
 * @return
 */
template<typename T>
[[nodiscard]] inline bool OrientationClockwise(const Polygon<T> &polygon) {
    return Area2(polygon) < 0.0;
}

/**
 * The same as clip::Orientation(clip::Geom2Clipper(polygon)).
 *
 * @tparam T
 * @param polygon
 * @return
 */
template<typename T>
[[nodiscard]] inline bool OrientationCounterClockwise(const Polygon<T> &polygon) {
    return Area2(polygon) >= 0.0;
}

template<typename T>
void ChangeOrientation(Polygon<T> &polygon) {
    std::reverse(polygon.begin(), polygon.end());
}

template<typename T>
[[nodiscard]] Polygon<T> ChangeOrientation(const Polygon<T> &polygon) {
    Polygon<T> ret(polygon.begin(), polygon.end());
    ChangeOrientation(ret);
    return ret;
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
std::vector<Tout> SampleSegment1D(
    const Tin &p0,
    const Tin &p1,
    const Tdist &max_sample_dist,
    Tdist &real_sample_dist,
    bool include_p0 = true,
    bool include_p1 = true
) {
    std::vector<Tout> sample_points;
    auto d = std::abs(static_cast<Tdist>(p0 - p1));
    int n_samp = (d == Tdist(0)) ? 2 : 1 + static_cast<int>(std::ceil(d / max_sample_dist));
    real_sample_dist = d / static_cast<Tdist>(n_samp - 1);
    sample_points.reserve(((include_p0) ? 1 : 0) + (n_samp - 2) + ((include_p1) ? 1 : 0));
    if (include_p0) {
        sample_points.push_back(static_cast<Tout>(p0));
    }
    if (n_samp > 2) {
        auto new_sample = static_cast<Tdist>(p0);
        for (int i = 0; i < n_samp - 2; ++i) {
            new_sample = new_sample + real_sample_dist;
            sample_points.push_back(static_cast<Tout>(new_sample));
        }
    }
    if (include_p1) {
        sample_points.push_back(static_cast<Tout>(p1));
    }
    return sample_points;
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
[[nodiscard]] inline std::vector<Tout> SampleSegment1D(
    const Tin &p0,
    const Tin &p1,
    const Tdist &max_sample_dist,
    bool include_p0 = true,
    bool include_p1 = true
) {
    Tdist unused;
    return SampleSegment1D<Tout>(p0, p1, max_sample_dist, unused, include_p0, include_p1);
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
Points<Tout> SampleSegment2D(
    const Point<Tin> &p0,
    const Point<Tin> &p1,
    const Tdist &max_sample_dist,
    Tdist &real_sample_dist,
    bool include_p0 = true,
    bool include_p1 = true
) {
    Points<Tin> sample_points;
    auto d = p0.template DistanceTo<Tdist>(p1);
    int n_samp = (d == Tdist(0)) ? 2 : 1 + static_cast<int>(std::ceil(d / max_sample_dist));
    real_sample_dist = d / static_cast<Tdist>(n_samp - 1);
    sample_points.reserve(((include_p0) ? 1 : 0) + (n_samp - 2) + ((include_p1) ? 1 : 0));
    if (include_p0) {
        sample_points.emplace_back(static_cast<Tout>(p0.x), static_cast<Tout>(p0.y));
    }
    if (n_samp > 2) {
        auto u = MakePoint(static_cast<Tdist>(p1.x - p0.x), static_cast<Tdist>(p1.y - p0.y));
        u = (u / u.Norm()) * real_sample_dist;
        auto new_sample = MakePoint(static_cast<Tdist>(p0.x), static_cast<Tdist>(p0.y));
        for (int i = 0; i < n_samp - 2; ++i) {
            new_sample = new_sample + u;
            sample_points.emplace_back(static_cast<Tout>(new_sample.x), static_cast<Tout>(new_sample.y));
        }
    }
    if (include_p1) {
        sample_points.emplace_back(static_cast<Tout>(p1.x), static_cast<Tout>(p1.y));
    }
    return sample_points;
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
[[nodiscard]] inline Points<Tout> SampleSegment2D(
    const Point<Tin> &p0,
    const Point<Tin> &p1,
    Tdist max_sample_dist,
    bool include_p0 = true,
    bool include_p1 = true
) {
    Tdist unused;
    return SampleSegment2D<Tout>(p0, p1, max_sample_dist, unused, include_p0, include_p1);
}

template<typename Tin, typename Tdist>
Tdist PolylineLength(
    const Points<Tin> &polyline,
    std::vector<Tdist> &accum_lengths
) {
    Tdist length = Tdist(0);
    accum_lengths.push_back(length);
    for (int i = 1; i < polyline.size(); ++i) {
        length += polyline[i - 1].template DistanceTo<Tdist>(polyline[i]);
        accum_lengths.push_back(length);
    }
    return length;
}

template<typename Tin, typename Tdist>
[[nodiscard]] inline Tdist PolylineLength(
    const Points<Tin> &polyline
) {
    std::vector<Tdist> unused;
    return PolylineLength(polyline, unused);
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
Points<Tout> SamplePolyline(
    const Points<Tin> &polyline,
    Tdist max_sample_dist,
    Tdist &real_sample_dist,
    bool include_front = true,
    bool include_back = true
) {
    std::vector<Tdist> accum_lengths;
    auto start = Tdist(0), end = PolylineLength(polyline, accum_lengths);
    auto sample_lengths = SampleSegment1D<Tdist>(start, end, max_sample_dist, real_sample_dist, false, false);
    int n_samples = static_cast<int>(sample_lengths.size());
    Points<Tout> sample_points;
    auto p_curr = geom::MakePoint(static_cast<Tout>(polyline.front().x), static_cast<Tout>(polyline.front().y));
    if (include_front) {
        sample_points.push_back(p_curr);
    }
    int sample_length_curr_idx = 0;
    for (int i = 1; i < polyline.size(); ++i) {
        const auto &p0 = polyline[i - 1];
        const auto &p1 = polyline[i];
        auto u = MakePoint(static_cast<Tdist>(p1.x - p0.x), static_cast<Tdist>(p1.y - p0.y));
        u = (u / u.Norm());
        while (sample_length_curr_idx < n_samples && accum_lengths[i] >= sample_lengths[sample_length_curr_idx]) {
            auto p = u * (sample_lengths[sample_length_curr_idx] - accum_lengths[i - 1]);
            sample_points.emplace_back(static_cast<Tout>(p0.x + p.x), static_cast<Tout>(p0.y + p.y));
            ++sample_length_curr_idx;
        }
    }
    if (include_back) {
        sample_points.emplace_back(static_cast<Tout>(polyline.back().x), static_cast<Tout>(polyline.back().y));
    }
    return sample_points;
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
[[nodiscard]] inline Points<Tout> SamplePolyline(
    const Points<Tin> &polyline,
    Tdist max_sample_dist,
    bool include_front = true,
    bool include_back = true
) {
    Tdist unused;
    return SamplePolyline<Tout>(polyline, max_sample_dist, unused, include_front, include_back);
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
inline Points<Tout> SamplePolygon(
    const Polygon<Tin> &polygon,
    Tdist max_sample_dist,
    Tdist &real_sample_dist
) {
    auto polygon_polyline = polygon;
    polygon_polyline.push_back(polygon.front());
    return SamplePolyline(polygon_polyline, max_sample_dist, real_sample_dist, true, false);
}

template<typename Tin, typename Tdist, typename Tout = Tdist>
[[nodiscard]] inline Points<Tout> SamplePolygon(
    const Polygon<Tin> &polygon,
    Tdist max_sample_dist
) {
    Tdist unused;
    return SamplePolygon<Tout>(polygon, max_sample_dist, unused);
}

template<typename T>
int Sgn(const T &val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T>
std::vector<T> SampleAngle(
    T a0,
    T a1,
    T max_sample_angle_abs,
    T &real_sample_angle,
    bool include_a0,
    bool include_a1
) {
    static_assert(std::is_floating_point<T>::value);
    std::vector<T> sample_angles;
    auto angle_diff = a1 - a0;
    auto angle_diff_abs = std::abs(angle_diff);
    int n_samp;
    if (angle_diff_abs < T(M_PIl)) {
        n_samp = (angle_diff == T(0)) ? 2 : static_cast<int>(std::ceil(angle_diff_abs / max_sample_angle_abs)) + 1;
        real_sample_angle = angle_diff / static_cast<T>(n_samp - 1);
    } else {
        angle_diff_abs = 2 * T(M_PIl) - angle_diff_abs;
        angle_diff = Sgn(angle_diff) * angle_diff_abs;
        n_samp = (angle_diff == T(0)) ? 2 : static_cast<int>(std::ceil(angle_diff_abs / max_sample_angle_abs)) + 1;
        real_sample_angle = -angle_diff / static_cast<T>(n_samp - 1);
    }
    sample_angles.reserve(((include_a0) ? 1 : 0) + (n_samp - 2) + ((include_a1) ? 1 : 0));
    if (include_a0) {
        sample_angles.push_back(a0);
    }
    if (n_samp > 2) {
        for (int i = 0; i < n_samp - 2; ++i) {
            sample_angles.push_back(((!i) ? a0 : sample_angles.back()) + real_sample_angle);
        }
    }
    if (include_a1) {
        sample_angles.push_back(a1);
    }
    return sample_angles;
}

template<typename T>
inline std::vector<T> SampleAngle(
    T a0,
    T a1,
    T max_sample_angle_abs,
    bool include_a0,
    bool include_a1
) {
    T unused;
    return SampleAngle(a0, a1, max_sample_angle_abs, unused, include_a0, include_a1);
}

template<typename T>
T MinimalTurningAngleAbs(
    const Point<T> &p_prev,
    const Point<T> &p,
    const Point<T> &p_next
) {
    auto u = p - p_prev;
    auto v = p_next - p;
    auto angle = std::abs(std::atan2(v.y, v.x) - std::atan2(u.y, u.x));
    angle = std::min(angle, T(2.0) * T(M_PIl) - angle);
    return angle;
}

template<typename T>
Point<T> ComputeConvexPolygonCenter(
    const Polygon<T> &polygon
) {
    auto center = MakePoint(T(0), T(0));
    double length_total = 0, length;
    for (int i_prev = static_cast<int>(polygon.size()) - 1, i = 0; i < polygon.size(); i_prev = i++) {
        const auto &p1 = polygon[i_prev];
        const auto &p2 = polygon[i];
        length = p1.DistanceTo(p2);
        length_total += length;
        center = center + (p1 + p2) * (length / 2.0);
    }
    center = center / length_total;
    return center;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Polygon<T> &polygon) {
    os << "{\n";
    for (unsigned j = 0; j < polygon.size(); ++j) {
        os << "  " << polygon[j] << (j + 1 == polygon.size() ? "" : ",") << "\n";
    }
    os << "}";
    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Polygons<T> &polygons) {
    os << "{\n";
    for (unsigned i = 0; i < polygons.size(); ++i) {
        os << "  {\n";
        for (unsigned j = 0; j < polygons[i].size(); ++j) {
            os << "    " << polygons[i][j] << (j + 1 == polygons[i].size() ? "" : ",") << "\n";
        }
        os << "  }" << (i + 1 == polygons.size() ? "" : ",") << "\n";
    }
    os << "}";
    return os;
}

}

#endif //TRIVIS_GENERIC_GEOM_UTILS_H_
