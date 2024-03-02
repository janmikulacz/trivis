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

#include "trivis/geom/poly_map.h"
#include "trivis/visibility_regions.h"

#include "trivis/mesh/tri_mesh.h"

#include "trivis/pl/point_location.h"

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

    struct ParamDefault {
        static constexpr double pl_bucket_size = 1.0;
        static constexpr bool pl_optimize_bucket_triangles = false;
        static constexpr std::array<double, 10> pl_eps1 = {1e-18, 1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9};
        static constexpr double pl_eps2_squared = 1e-12 * 1e-12;
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
        int triangle_id;
        std::vector<int> snap_to_nodes;
    };

    struct RayShootingResult {
        char code = 'X';
        geom::FPoint p;
        geom::FPoint p2;
        std::optional<int> edge_id;
        std::optional<int> node_id;
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

    std::optional<bool> IsVisible(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<bool> IsVisible(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<bool> IsVisible(
        int node_id,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: INTERSECTION OF RAY AND OBSTACLE ##### ///

    std::optional<RayShootingResult> ShootRay(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const geom::FPoint &direction,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<RayShootingResult> ShootRay(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoint &direction,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<RayShootingResult> ShootRay(
        int node_id,
        const geom::FPoint &direction,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: MAP VERTICES ##### ///

    std::optional<std::vector<int>> VisibleVertices(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> VisibleVertices(
        const geom::FPoint &q,
        int q_triangle_id,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> VisibleVertices(
        int node_id,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: INPUT POINTS ##### ///

    std::optional<std::vector<int>> VisiblePoints(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> VisiblePoints(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> VisiblePoints(
        int node_id,
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> VisiblePoints(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> VisiblePoints(
        int node_id,
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: VISIBILITY GRAPHS ##### ///

    std::optional<std::vector<std::vector<int>>> VertexVertexVisibilityGraph(
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<bool>>> VertexVertexVisibilityGraphBool(
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<int>>> VertexPointVisibilityGraph(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<bool>>> VertexPointVisibilityGraphBool(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<int>>> PointPointVisibilityGraph(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<bool>>> PointPointVisibilityGraphBool(
        const geom::FPoints &points,
        const std::vector<std::optional<PointLocationResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: VISIBILITY REGIONS (ABSTRACT REPRESENTATION) ##### ///

    std::optional<AbstractVisibilityRegion> VisibilityRegion(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegion(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegion(
        int node_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegionIterative(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegionIterative(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegionIterative(
        int node_id,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegionWithHistory(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::vector<ExpansionHistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegionWithHistory(
        const geom::FPoint &q,
        int q_triangle_id,
        std::vector<ExpansionHistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> VisibilityRegionWithHistory(
        int node_id,
        std::vector<ExpansionHistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### VISIBILITY REGIONS POSTPROCESSING AND UTILITIES ##### ///

    [[nodiscard]] RadialVisibilityRegion ComputeIntersections(
        const AbstractVisibilityRegion &abstract,
        bool line_line_mode = true
    ) const;

    [[nodiscard]] static bool IsValid(const RadialVisibilityRegion &visibility_region);

    static void RemoveAntennas(
        RadialVisibilityRegion &visibility_region
    );

    [[nodiscard]] static RadialVisibilityRegion IntersectWithCircle(
        double radius,
        const RadialVisibilityRegion &visibility_region
    );

    static void RemoveShortEdges(
        double min_edge_length,
        RadialVisibilityRegion &visibility_region
    );

    [[nodiscard]] static RadialVisibilityRegion SampleArcEdges(
        const RadialVisibilityRegion &visibility_region,
        double max_angle
    );

    [[nodiscard]] RadialVisibilityRegion Postprocess(
        const AbstractVisibilityRegion &abstract,
        bool line_line_mode_intersections = true,
        bool remove_antennas = false,
        std::optional<double> radius = std::nullopt,
        std::optional<double> min_edge_length = std::nullopt,
        std::optional<double> sampling_max_angle = std::nullopt
    ) const;

    [[nodiscard]] static geom::FPolygon ConvertToPolygon(
        const RadialVisibilityRegion &visibility_region
    );

    /// ##### SHORTCUTS: VISIBILITY REGIONS + POSTPROCESSING ##### ///

    std::optional<RadialVisibilityRegion> VisibilityRegionWithPostprocessing(
        const geom::FPoint &q,
        const PointLocationResult &q_location,
        std::optional<double> radius = std::nullopt,
        bool line_line_mode_intersections = true,
        bool remove_antennas = false,
        std::optional<double> min_edge_length = std::nullopt,
        std::optional<double> sampling_max_angle = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<RadialVisibilityRegion> VisibilityRegionWithPostprocessing(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        bool line_line_mode_intersections = true,
        bool remove_antennas = false,
        std::optional<double> min_edge_length = std::nullopt,
        std::optional<double> sampling_max_angle = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    std::optional<RadialVisibilityRegion> VisibilityRegionWithPostprocessing(
        int node_id,
        std::optional<double> radius = std::nullopt,
        bool line_line_mode_intersections = true,
        bool remove_antennas = false,
        std::optional<double> min_edge_length = std::nullopt,
        std::optional<double> sampling_max_angle = std::nullopt,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### GENERAL UTILITIES ##### ///

    static bool SaveMesh(
        const mesh::TriMesh &mesh,
        const std::string &file,
        std::stringstream *info = nullptr
    );

    static std::optional<mesh::TriMesh> LoadMesh(
        const std::string &file,
        std::stringstream *info = nullptr
    );

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
    std::vector<double> _pl_eps1_seq{ParamDefault::pl_eps1.begin(), ParamDefault::pl_eps1.end()};
    double _pl_eps2_squared = ParamDefault::pl_eps2_squared;

    /// ##### EXPAND EDGE: INTERSECTION OF RAY AND OBSTACLE ##### ///

    [[nodiscard]] RayShootingResult ExpandEdgeObstacleIntersection(
        int level,
        const geom::FPoint &q,
        const geom::FPoint &distant_t,
        int node_l_id,
        int node_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        ExpansionStats *stats = nullptr
    ) const;

    /// ##### EXPAND EDGE: TWO-POINT QUERIES ##### ///

    [[nodiscard]] bool ExpandEdgeVisibilityBetween(
        int level,
        const geom::FPoint &q,
        const geom::FPoint &t,
        int node_l_id,
        int node_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
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
        const std::vector<bool> *tabu_vertices = nullptr,
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
