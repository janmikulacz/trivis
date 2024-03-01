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

#include "trivis/bucketing/bucketing.h"

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

    struct Stats {
        int num_expansions = 0;
        int max_recursion_depth = 0;
    };

    struct HistoryStep {
        int edge_id = -1;
        std::optional<int> tri_id = std::nullopt;
        int rest_l_id = -1;
        int rest_r_id = -1;
        std::optional<int> res_seg_id = std::nullopt;
    };

    struct ObstacleIntersection {
        char code = 'X';
        geom::FPoint p;
        geom::FPoint p2;
        std::optional<int> edge_id;
        std::optional<int> node_id;
    };

    struct LocatePointResult {
        int triangle_id;
        std::vector<int> snap_nodes;
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

    /// ##### GETTERS ##### ///

    [[nodiscard]] bool has_map() const {
        return _has_map;
    }

    [[nodiscard]] bool has_mesh() const {
        return _has_mesh;
    }

    [[nodiscard]] bool has_bucketing() const {
        return _has_mesh;
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

    [[nodiscard]] const auto &bucketing() const {
        assert(_has_bucketing);
        return _bucketing;
    }

    /// ##### INITIALIZATION ##### ///

    void SetMap(geom::PolyMap map);

    void ConstructMeshConstrainedDelaunayTriangulation();

    void SetMesh(
        mesh::TriMesh mesh,
        std::optional<geom::PolyMap> map = std::nullopt
    );

    void FillBucketing(
        std::optional<double> bucket_size = std::nullopt,
        std::optional<double> mean_bucket_triangle_count_max = std::nullopt
    );

    void OptimizeBuckets();

    void SetBucketingExact(std::optional<bool> exact = std::nullopt) {
        assert(_has_bucketing);
        _bucketing.SetExact(exact);
    }

    void SetBucketingEpsilons(const std::optional<std::vector<double>> &epsilons = std::nullopt) {
        assert(_has_bucketing);
        _bucketing.SetEpsilons(epsilons);
    }

    /// ##### LOCATE POINT ##### ///

    [[nodiscard]] std::optional<LocatePointResult> LocatePoint(
        const geom::FPoint &q,
        double epsilon = 1e-24
    ) const;

    /// ##### VISIBILITY: INTERSECTION OF RAY AND OBSTACLE ##### ///

    std::optional<ObstacleIntersection> ComputeObstacleIntersection(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        int q_triangle_id,
        const geom::FPoint &direction,
        Stats *stats = nullptr
    ) const;

    std::optional<ObstacleIntersection> ComputeObstacleIntersection(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoint &direction,
        Stats *stats = nullptr
    ) const;

    std::optional<ObstacleIntersection> ComputeObstacleIntersection(
        int node_id,
        const geom::FPoint &direction,
        Stats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: TWO-POINT QUERIES ##### ///

    std::optional<bool> ComputeVisibilityBetween(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<bool> ComputeVisibilityBetween(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<bool> ComputeVisibilityBetween(
        int node_id,
        const geom::FPoint &p,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: MAP VERTICES ##### ///

    std::optional<std::vector<int>> ComputeVisibleVertices(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> ComputeVisibleVertices(
        const geom::FPoint &q,
        int q_triangle_id,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> ComputeVisibleVertices(
        int node_id,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: INPUT POINTS ##### ///

    std::optional<std::vector<int>> ComputeVisiblePoints(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> ComputeVisiblePoints(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> ComputeVisiblePoints(
        int node_id,
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> ComputeVisiblePoints(
        const geom::FPoint &q,
        int q_triangle_id,
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<int>> ComputeVisiblePoints(
        int node_id,
        const geom::FPoints &points,
        const std::vector<std::vector<int>> &triangle_points,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: VISIBILITY GRAPHS ##### ///

    std::optional<std::vector<std::vector<int>>> ComputeVisibilityGraphNodeNode(
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<bool>>> ComputeVisibilityGraphNodeNodeBoolMatrix(
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<int>>> ComputeVisibilityGraphPointNode(
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<bool>>> ComputeVisibilityGraphPointNodeBoolMatrix(
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        const std::vector<bool> *tabu_vertices = nullptr,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<int>>> ComputeVisibilityGraphPointPoint(
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<std::vector<std::vector<bool>>> ComputeVisibilityGraphPointPointBoolMatrix(
        const geom::FPoints &points,
        const std::vector<std::optional<LocatePointResult>> &points_locations,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    /// ##### VISIBILITY: VISIBILITY REGIONS (ABSTRACT REPRESENTATION) ##### ///

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegion(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegion(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegion(
        int node_id,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegionIterative(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegionIterative(
        const geom::FPoint &q,
        int q_triangle_id,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegionIterative(
        int node_id,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegionWithHistory(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        std::vector<HistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegionWithHistory(
        const geom::FPoint &q,
        int q_triangle_id,
        std::vector<HistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    std::optional<AbstractVisibilityRegion> ComputeVisibilityRegionWithHistory(
        int node_id,
        std::vector<HistoryStep> &history,
        std::optional<double> radius = std::nullopt,
        Stats *stats = nullptr
    ) const;

    /// ##### VISIBILITY REGIONS POSTPROCESSING AND UTILITIES ##### ///

    [[nodiscard]] RadialVisibilityRegion ComputeIntersections(
        const AbstractVisibilityRegion &abstract,
        bool fast_mode = true
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
        double max_sample_beta
    );

    [[nodiscard]] static geom::FPolygon ConvertToPolygon(
        const RadialVisibilityRegion &visibility_region
    );

    /// ##### SHORTCUTS: VISIBILITY REGIONS + POSTPROCESSING ##### ///

    std::optional<RadialVisibilityRegion> ComputeRadialVisibilityRegion(
        const geom::FPoint &q,
        const LocatePointResult &q_location,
        std::optional<double> radius,
        bool intersections_fast_mode = true,
        bool remove_antennas = true,
        std::optional<double> min_edge_length = std::nullopt,
        Stats *stats = nullptr
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
    bool _has_bucketing = false;
    geom::FLimits _limits;
    geom::PolyMap _map;
    mesh::TriMesh _mesh;
    geom::FPolygons _triangles;
    bucketing::Bucketing _bucketing;

    /// ##### EXPAND EDGE: INTERSECTION OF RAY AND OBSTACLE ##### ///

    [[nodiscard]] ObstacleIntersection ExpandEdgeObstacleIntersection(
        int level,
        const geom::FPoint &q,
        const geom::FPoint &distant_t,
        int node_l_id,
        int node_r_id,
        int curr_edge_id,
        int curr_edge_tri_id,
        Stats *stats = nullptr
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
        Stats *stats = nullptr
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
        Stats *stats = nullptr
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
        Stats *stats = nullptr
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
        Stats *stats = nullptr
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
        Stats *stats = nullptr
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
        std::vector<HistoryStep> &history,
        Stats *stats = nullptr
    ) const;

    /// ##### POSTPROCESSING HELPER METHODS ##### ///

    void AppendNotCollinearWithIntersection(
        const geom::FPoint &q,
        const AbstractVisibilityRegionVertex &v,
        int edge_flag,
        bool is_last,
        bool fast_mode,
        RadialVisibilityRegion &visibility_region
    ) const;

};

}

#endif //TRIVIS_TRIVIS_H_
