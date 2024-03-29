/**
 * File:   trivis.h
 *
 * Date:   29.11.2021
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#ifndef TRIVIS_TRIVIS_H_
#define TRIVIS_TRIVIS_H_

#include "trivis/geom/generic_geom_types.h"
#include "trivis/geom/generic_geom_utils.h"
#include "trivis/geom/geom_types.h"
#include "trivis/geom/poly_map.h"
#include "trivis/geom/robust_geometry.h"

#include "trivis/mesh/cdt.h"
#include "trivis/mesh/tri_mesh.h"

#include "trivis/pl/grid.h"
#include "trivis/pl/point_location.h"

#include "trivis/utils/clipper_geom.h"
#include "trivis/utils/clipper_utils.h"
#include "trivis/utils/generic_utils.h"
#include "trivis/utils/random_points.h"
#include "trivis/utils/simple_clock.h"

#include "trivis/vis_regions.h"

#include <memory>
#include <utility>
#include <cassert>

namespace trivis {

class Trivis {

public:

    ///================///
    /// PUBLIC MEMBERS ///
    ///================///

    /// ##### STRUCTURES ##### ///

    struct DefaultParam {
        static constexpr double kPointLocationBucketSize = 1.0;
        static constexpr bool kPointLocationOptimizeBucketTriangles = false;
        static constexpr std::array<double, 10> kPointLocationEpsilon1 = {1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9};
        static constexpr double kPointLocationEpsilon2Squared = 1e-12 * 1e-12;
    };

    struct ExpansionStats {
        int num_expansions = 0;
        int max_recursion_depth = 0;
    };

    struct ExpansionHistoryStep {
        int edge_id = -1;
        int rest_l_id = -1;
        int rest_r_id = -1;
        std::optional<int> tri_id = std::nullopt;
        std::optional<int> output_segment_id = std::nullopt;
    };

    struct PointLocationResult {
        int tri_id;
        std::vector<int> snap_to_vertices;
    };

    struct RayShootingResult {
        char code = 'X';
        geom::FPoint p;
        geom::FPoint p2;
        std::optional<int> edge_id;
        std::optional<int> ver_id;
    };


    /// ##### CONSTRUCTORS ##### ///

    /// Default constructor.
    Trivis() = default;

    /// Copy constructor.
    Trivis(const Trivis &) = default;

    /// Move constructor.
    Trivis(Trivis &&) noexcept = default;

    /// Copy operator.
    Trivis &operator=(const Trivis &) = default;

    /// Move operator.
    Trivis &operator=(Trivis &&) noexcept = default;

    /// Destructor.
    virtual ~Trivis() = default;

    /// Constructor with map initialization.
    explicit Trivis(
        geom::PolyMap map,
        std::optional<double> pl_bucket_size = std::nullopt,
        std::optional<double> pl_max_avg_triangle_count_in_bucket = std::nullopt,
        std::optional<bool> pl_optimize_bucket_triangles = std::nullopt,
        const std::optional<std::vector<double>> &pl_eps1_seq = std::nullopt,
        std::optional<double> pl_eps2_squared = std::nullopt
    ) {
        Init(std::move(map), pl_bucket_size, pl_max_avg_triangle_count_in_bucket, pl_optimize_bucket_triangles, pl_eps1_seq, pl_eps2_squared);
    }

    /// Constructor with mesh initialization.
    explicit Trivis(
        mesh::TriMesh mesh,
        std::optional<geom::PolyMap> map = std::nullopt,
        std::optional<double> pl_bucket_size = std::nullopt,
        std::optional<double> pl_max_avg_triangle_count_in_bucket = std::nullopt,
        std::optional<bool> pl_optimize_bucket_triangles = std::nullopt,
        const std::optional<std::vector<double>> &pl_eps1_seq = std::nullopt,
        std::optional<double> pl_eps2_squared = std::nullopt
    ) {
        Init(std::move(mesh), std::move(map), pl_bucket_size, pl_max_avg_triangle_count_in_bucket, pl_optimize_bucket_triangles, pl_eps1_seq, pl_eps2_squared);
    }

    /// ##### GETTERS ##### ///

    [[nodiscard]] bool has_map() const {
        return _has_map;
    }

    [[nodiscard]] bool has_mesh() const {
        return _has_mesh;
    }

    [[nodiscard]] bool has_pl() const {
        return _has_pl;
    }

    [[nodiscard]] const auto &limits() const {
        assert(_has_map || _has_mesh);
        return _limits;
    }

    [[nodiscard]] const auto &map() const {
        assert(_has_map);
        return _map;
    }

    [[nodiscard]] const auto &mesh() const {
        assert(_has_mesh);
        return _mesh;
    }

    [[nodiscard]] const auto &triangles() const {
        assert(_has_mesh);
        return _triangles;
    }

    [[nodiscard]] const auto &pl() const {
        assert(_has_pl);
        return _pl;
    }

    /// ##### INITIALIZATION ##### ///

    void SetMap(geom::PolyMap map);

    void ConstructMeshCDT();

    void SetMesh(
        mesh::TriMesh mesh,
        std::optional<geom::PolyMap> map = std::nullopt
    );

    void FillPointLocationBuckets(
        std::optional<double> bucket_size = std::nullopt,
        std::optional<double> max_avg_triangle_count_in_bucket = std::nullopt
    );

    void OptimizePointLocationBucketTriangles();

    void InitPointLocation(
        std::optional<double> bucket_size,
        std::optional<double> max_avg_triangle_count_in_bucket = std::nullopt,
        std::optional<bool> optimize_bucket_triangles = std::nullopt
    );

    void SetPointLocationEpsilons(
        const std::optional<std::vector<double>> &eps1_seq = std::nullopt,
        std::optional<double> eps2_squared = std::nullopt
    );

    void Init(
        geom::PolyMap map,
        std::optional<double> pl_bucket_size = std::nullopt,
        std::optional<double> pl_max_avg_triangle_count_in_bucket = std::nullopt,
        std::optional<bool> pl_optimize_bucket_triangles = std::nullopt,
        const std::optional<std::vector<double>> &pl_eps1_seq = std::nullopt,
        std::optional<double> pl_eps2_squared = std::nullopt
    );

    void Init(
        mesh::TriMesh mesh,
        std::optional<geom::PolyMap> map = std::nullopt,
        std::optional<double> pl_bucket_size = std::nullopt,
        std::optional<double> pl_max_avg_triangle_count_in_bucket = std::nullopt,
        std::optional<bool> pl_optimize_bucket_triangles = std::nullopt,
        const std::optional<std::vector<double>> &pl_eps1_seq = std::nullopt,
        std::optional<double> pl_eps2_squared = std::nullopt
    );

    /// ##### LOCATE POINT ##### ///

    [[nodiscard]] std::optional<PointLocationResult> LocatePoint(
        const geom::FPoint &q,
        const std::optional<std::vector<double>> &eps1_seq = std::nullopt,
        const std::optional<double> &eps2_squared = std::nullopt
    ) const;

    /// ##### VISIBILITY: TWO-POINT QUERIES ##### ///

    bool IsVisible(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    bool IsVisible(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    bool IsVisible(
        int ver_id,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: INTERSECTION OF RAY AND OBSTACLE ##### ///

    std::optional<RayShootingResult> ShootRay(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const geom::FPoint &direction,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<RayShootingResult> ShootRay(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoint &direction,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<RayShootingResult> ShootRay(
        int ver_id,
        const geom::FPoint &direction,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: MAP VERTICES ##### ///

    std::vector<int> VisibleVertices(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<int> VisibleVertices(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<int> VisibleVertices(
        int ver_id,
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: INPUT POINTS ##### ///

    std::vector<int> VisiblePoints(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<int> VisiblePoints(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<int> VisiblePoints(
        int ver_id,
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<int> VisiblePoints(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<int> VisiblePoints(
        int ver_id,
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: VISIBILITY GRAPHS ##### ///

    std::vector<std::vector<int>> VertexVertexVisibilityGraph(
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<std::vector<bool>> VertexVertexVisibilityGraphBool(
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<std::vector<int>> PointVertexVisibilityGraph(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<std::vector<bool>> PointVertexVisibilityGraphBool(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<std::vector<int>> PointPointVisibilityGraph(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::vector<std::vector<bool>> PointPointVisibilityGraphBool(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: VISIBILITY REGIONS (ABSTRACT REPRESENTATION) ##### ///

    AbstractVisibilityRegion VisibilityRegion(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegion(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegion(
        int ver_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegionIterative(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegionIterative(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegionIterative(
        int ver_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegionWithHistory(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::vector<ExpansionHistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegionWithHistory(
        const geom::FPoint &q,
        int q_triangle_id,
        std::vector<ExpansionHistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    AbstractVisibilityRegion VisibilityRegionWithHistory(
        int ver_id,
        std::vector<ExpansionHistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY REGIONS POSTPROCESSING AND UTILITIES ##### ///

    [[nodiscard]] RadialVisibilityRegion ToRadialVisibilityRegion(
        const AbstractVisibilityRegion &abstract,
        bool line_line_mode = true
    ) const;

private:

    ///=================///
    /// PRIVATE MEMBERS ///
    ///=================///

    /// ##### STRUCTURES ##### ///

    struct EdgeExpansionInfo {
        int level = -1;
        int rest_l_id = -1;
        int rest_r_id = -1;
        int curr_edge_id = -1;
        int curr_edge_tri_id = -1;
    };

    /// ##### MEMBER VARIABLES ##### ///

    bool _has_map = false;
    bool _has_mesh = false;
    bool _has_pl = false;
    geom::FLimits _limits;
    geom::PolyMap _map;
    mesh::TriMesh _mesh;
    geom::FPolygons _triangles;
    pl::PointLocation _pl;
    std::vector<double> _pl_eps1_seq{DefaultParam::kPointLocationEpsilon1.begin(), DefaultParam::kPointLocationEpsilon1.end()};
    double _pl_eps2_squared = DefaultParam::kPointLocationEpsilon2Squared;

    /// ##### EXPAND EDGE: TWO-POINT QUERIES ##### ///

    [[nodiscard]] bool ExpandEdge2PointVisibility(
        int level,
        const geom::FPoint &q,
        const geom::FPoint &t,
        int rest_l_id,
        int rest_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### EXPAND EDGE: INTERSECTION OF RAY AND OBSTACLE ##### ///

    [[nodiscard]] std::optional<RayShootingResult> ExpandEdgeRayShooting(
        int level,
        const geom::FPoint &q,
        const geom::FPoint &distant_t,
        int rest_l_id,
        int rest_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        double sq_radius,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### EXPAND EDGE: MAP VERTICES ##### ///

    void ExpandEdgeVisibleVertices(
        int level,
        const geom::FPoint &q,
        int rest_l_id,
        int rest_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        double sq_radius,
        std::vector<int> &visible_vertices,
        const std::optional<std::vector<bool>> &tabu_vertices = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### EXPAND EDGE: INPUT POINTS ##### ///

    void ExpandEdgeVisiblePoints(
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        int level,
        const geom::FPoint &q,
        int rest_l_id,
        int rest_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        double sq_radius,
        std::vector<bool> &point_visited,
        std::vector<int> &visible_points,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### EXPAND EDGE: VISIBILITY REGIONS ##### ///

    void ExpandEdgeVisibilityRegion(
        int level,
        const geom::FPoint &q,
        int rest_l_id,
        int rest_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        double sq_radius,
        AbstractVisibilityRegion &visibility_region,
        ExpansionStats *stats = nullptr
    ) const;

    void ExpandEdgeVisibilityRegionIterative(
        const geom::FPoint &q,
        const EdgeExpansionInfo &expansion,
        double sq_radius,
        AbstractVisibilityRegion &visibility_region,
        bool &has_expansion1,
        EdgeExpansionInfo &expansion1,
        bool &has_expansion2,
        EdgeExpansionInfo &expansion2,
        ExpansionStats *stats = nullptr
    ) const;

    void ExpandEdgeVisibilityRegionWithHistory(
        int level,
        const geom::FPoint &q,
        int rest_l_id,
        int rest_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        double sq_radius,
        AbstractVisibilityRegion &visibility_region,
        std::vector<ExpansionHistoryStep> &history,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### POSTPROCESSING HELPER METHODS ##### ///

    void AppendNotCollinearWithIntersection(
        const geom::FPoint &q,
        const AbstractVisibilityRegionVertex &v,
        int edge_flag,
        bool is_last,
        bool line_line_mode,
        RadialVisibilityRegion &visibility_region
    ) const;

};

}

#endif //TRIVIS_TRIVIS_H_
