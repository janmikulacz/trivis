/**
 * File:   point_location.h
 *
 * Date:   24.03.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_PL_POINT_LOCATION_H_
#define TRIVIS_PL_POINT_LOCATION_H_

#include <array>
#include <optional>

#include "trivis/pl/grid.h"

#include "trivis/geom/geom_types.h"

namespace trivis::pl {

void FillBuckets(
    Grid<std::vector<int>> &buckets,
    const geom::FLimits &lim,
    double x_scale,
    double y_scale,
    const geom::FPolygons &triangles
);

void SortTrianglesInBucketsByLargestIntersectionArea(
    Grid<std::vector<int>> &buckets,
    const geom::FLimits &lim,
    double x_scale,
    double y_scale,
    const geom::FPolygons &triangles
);

void RemoveDuplicateTrianglesInBuckets(
    Grid<std::vector<int>> &buckets
);

std::optional<int> FindTriangleContainingPoint(
    const geom::FPoint &q,
    const geom::FPolygons &triangles,
    const Grid<std::vector<int>> &buckets,
    const geom::FLimits &lim,
    double x_scale,
    double y_scale,
    const std::vector<double> &epsilons
);

class PointLocation {
public:

    void Init(
        geom::FLimits lim,
        const geom::FPolygons &triangles,
        double bucket_size
    ) {
        _lim = lim;
        _n_col = static_cast<int>(std::ceil((lim.x_max - lim.x_min) / bucket_size));
        _n_row = static_cast<int>(std::ceil((lim.y_max - lim.y_min) / bucket_size));
        _x_scale = (lim.x_max - lim.x_min) / static_cast<double>(_n_col);
        _y_scale = (lim.y_max - lim.y_min) / static_cast<double>(_n_row);
        _buckets.Clear();
        _buckets.Resize(_n_row, _n_col);
        pl::FillBuckets(_buckets, _lim, _x_scale, _y_scale, triangles);
    }

    void SortTrianglesInBucketsByLargestIntersectionArea(const geom::FPolygons &triangles) {
        pl::SortTrianglesInBucketsByLargestIntersectionArea(_buckets, _lim, _x_scale, _y_scale, triangles);
    }

    void RemoveDuplicateTrianglesInBuckets() {
        pl::RemoveDuplicateTrianglesInBuckets(_buckets);
    }

    [[nodiscard]] std::optional<int> FindTriangle(
        const geom::FPoint &q,
        const geom::FPolygons &triangles,
        const std::vector<double> &epsilons
    ) const {
        return pl::FindTriangleContainingPoint(q, triangles, _buckets, _lim, _x_scale, _y_scale, epsilons);
    }

    [[nodiscard]] int ToIdxX(double x) const;

    [[nodiscard]] int ToIdxY(double y) const;

    [[nodiscard]] const geom::FLimits &lim() const {
        return _lim;
    }

    [[nodiscard]] int n_row() const {
        return _n_row;
    }

    [[nodiscard]] int n_col() const {
        return _n_col;
    }

    [[nodiscard]] double x_scale() const {
        return _x_scale;
    }

    [[nodiscard]] double y_scale() const {
        return _y_scale;
    }

    [[nodiscard]] const Grid<std::vector<int>> &buckets() const {
        return _buckets;
    }

private:
    geom::FLimits _lim;
    int _n_row = 0;
    int _n_col = 0;
    double _x_scale = 0.0;
    double _y_scale = 0.0;
    Grid<std::vector<int>> _buckets;
};

}

#endif //TRIVIS_PL_POINT_LOCATION_H_
