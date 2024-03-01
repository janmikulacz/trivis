/**
 * File:   bucketing.h
 *
 * Date:   24.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_BUCKETING_BUCKETING_H_
#define TRIVIS_BUCKETING_BUCKETING_H_

#include <array>
#include <optional>

#include "trivis/bucketing/grid.h"

#include "trivis/geom/geom_types.h"

namespace trivis::bucketing {

void FillBuckets(
    Grid &buckets,
    const geom::FLimits &lim,
    double x_scale,
    double y_scale,
    const geom::FPolygons &triangles
);

void OrderBucketTrianglesByIntersectionArea(
    Grid &buckets,
    const geom::FLimits &lim,
    double x_scale,
    double y_scale,
    const geom::FPolygons &triangles
);

void RemoveDuplicateBucketTriangles(
    Grid &buckets
);

std::optional<int> FindTriangle(
    const geom::FPoint &q,
    const geom::FPolygons &triangles,
    const Grid &buckets,
    const geom::FLimits &lim,
    double x_scale,
    double y_scale,
    bool exact,
    const std::vector<double> &epsilons
);

class Bucketing {
public:

    static constexpr double kDefaultBucketSize = 1.0;
    static constexpr bool kDefaultExact = false;
    static constexpr std::array<double, 10> kDefaultEpsilons = {1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9};

    void Init(
        geom::FLimits lim,
        double bucket_size = kDefaultBucketSize
    ) {
        _lim = lim;
        _ncol = static_cast<int>(std::ceil((lim.x_max - lim.x_min) / bucket_size));
        _nrow = static_cast<int>(std::ceil((lim.y_max - lim.y_min) / bucket_size));
        _x_scale = (lim.x_max - lim.x_min) / static_cast<double>(_ncol);
        _y_scale = (lim.y_max - lim.y_min) / static_cast<double>(_nrow);
        _buckets.Resize(_nrow, _ncol);
    }

    void FillBuckets(const geom::FPolygons &triangles) {
        _buckets.Clear();
        bucketing::FillBuckets(_buckets, _lim, _x_scale, _y_scale, triangles);
    }

    void OrderBucketTrianglesByIntersectionArea(const geom::FPolygons &triangles) {
        bucketing::OrderBucketTrianglesByIntersectionArea(_buckets, _lim, _x_scale, _y_scale, triangles);
    }

    void RemoveDuplicateBucketTriangles() {
        bucketing::RemoveDuplicateBucketTriangles(_buckets);
    }

    void SetExact(std::optional<bool> exact) {
        _exact = exact.value_or(kDefaultExact);
    }

    void SetEpsilons(const std::optional<std::vector<double>> &epsilons) {
        _epsilons = epsilons.value_or(std::vector<double>(kDefaultEpsilons.begin(), kDefaultEpsilons.end()));
    }

    [[nodiscard]] std::optional<int> FindTriangle(
        const geom::FPoint &q,
        const geom::FPolygons &triangles,
        std::optional<bool> exact = std::nullopt
    ) const {
        return bucketing::FindTriangle(q, triangles, _buckets, _lim, _x_scale, _y_scale, exact.value_or(_exact), _epsilons);
    }

    [[nodiscard]] int ToIdxX(double x) const;

    [[nodiscard]] int ToIdxY(double y) const;

    [[nodiscard]] const geom::FLimits &lim() const {
        return _lim;
    }
    [[nodiscard]] int nrow() const {
        return _nrow;
    }
    [[nodiscard]] int ncol() const {
        return _ncol;
    }
    [[nodiscard]] double x_scale() const {
        return _x_scale;
    }
    [[nodiscard]] double y_scale() const {
        return _y_scale;
    }

    [[nodiscard]] const Grid &buckets() const {
        return _buckets;
    }

    [[nodiscard]] bool exact() const {
        return _exact;
    }

    [[nodiscard]] const std::vector<double> &epsilons() const {
        return _epsilons;
    }

private:
    geom::FLimits _lim;
    int _nrow = 0;
    int _ncol = 0;
    double _x_scale = 0.0;
    double _y_scale = 0.0;
    Grid _buckets;
    bool _exact = kDefaultExact;
    std::vector<double> _epsilons = std::vector<double>(kDefaultEpsilons.begin(), kDefaultEpsilons.end());
};

}

#endif //TRIVIS_BUCKETING_BUCKETING_H_
