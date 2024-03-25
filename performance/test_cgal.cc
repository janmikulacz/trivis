/**
 * File:   test_cgal.cc
 *
 * Date:   13.01.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "trivis/utils/random_points.h"
#include "trivis/utils/simple_clock.h"

#include "trivis_plus/data_loading/load_map.h"
#include "trivis_plus/utils/log.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>

// #include <CGAL/Arr_naive_point_location.h>
// #include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
// #include <CGAL/Arr_trapezoid_ric_point_location.h>
// #include <CGAL/Arr_triangulation_point_location.h>

#ifndef DEFAULT_MAP_DIR
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_POINT_DIR
#define DEFAULT_POINT_DIR "."
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;

template<typename Kernel>
using Point = typename Kernel::Point_2;

template<typename Kernel>
using Segment = typename Kernel::Segment_2;

template<typename Kernel>
using Traits = CGAL::Arr_segment_traits_2<Kernel>;

template<typename Kernel>
using Arrangement = CGAL::Arrangement_2<Traits<Kernel>>;

template<typename Kernel>
using Fh = typename Arrangement<Kernel>::Face_const_handle;

template<typename Kernel>
using Hh = typename Arrangement<Kernel>::Halfedge_const_handle;

template<typename Kernel>
using Vh = typename Arrangement<Kernel>::Vertex_const_handle;

template<typename Kernel>
using VisibilitySimplePolygon = CGAL::Simple_polygon_visibility_2<Arrangement<Kernel>, CGAL::Tag_false>; // CGAL::Tag_true to remove antennas

template<typename Kernel>
using VisibilityRotationalSweep = CGAL::Rotational_sweep_visibility_2<Arrangement<Kernel>, CGAL::Tag_false>; // CGAL::Tag_true to remove antennas

template<typename Kernel>
using VisibilityTriangularExpansion = CGAL::Triangular_expansion_visibility_2<Arrangement<Kernel>, CGAL::Tag_false>; // CGAL::Tag_true to remove antennas

// template<typename Kernel>
// using PointLocation = CGAL::Arr_naive_point_location<Arrangement<Kernel>>; // 9.47961 s (OK)

// template<typename Kernel>
// using PointLocation = CGAL::Arr_landmarks_point_location<Arrangement<Kernel>>; // 5.95374 s (OK)

template<typename Kernel>
using PointLocation = CGAL::Arr_walk_along_line_point_location<Arrangement<Kernel>>; // 2.83884 s (OK)

// template<typename Kernel>
// using PointLocation = CGAL::Arr_trapezoid_ric_point_location<Arrangement<Kernel>>; // 0.0474296 s (FAIL)

// template<typename Kernel>
// using PointLocation = CGAL::Arr_triangulation_point_location<Arrangement<Kernel>>; // 0.0399803 s (FAIL)

struct ProgramOptionVariables {
    std::string map_name;
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string points_name;
    std::string points_extension = ".txt";
    std::string points_dir = DEFAULT_POINT_DIR;
    std::string points_full_path;
    std::string out_file;
    bool resume_after_crash = false;
    bool check_map_validity = false;
    std::string algorithm = "triangular_expansion"; // "simple_polygon", "rotational_sweep", "triangular_expansion
    bool inexact_constructions = false;
    int max_points = std::numeric_limits<int>::max();
};

template<typename T>
inline auto po_value(T &variable) {
    return po::value(&variable)->default_value(variable);
}

void AddProgramOptions(
    po::options_description &options_description,
    ProgramOptionVariables &pov
) {
    options_description.add_options()
        ("help,h", "Produce this help message. \n (*) Overwrites options: all.")
        ("map-name", po_value(pov.map_name), "Map name.")
        ("map-extension", po_value(pov.map_extension), "Map extension.")
        ("map-dir", po_value(pov.map_dir), "Map directory.")
        ("map-full-path", po_value(pov.map_full_path), "Full path to the map file. \n (*) Overwrites options: map-name, map-ext, map-dir.")
        ("points-name", po_value(pov.points_name), "Points name.")
        ("points-extension", po_value(pov.points_extension), "Points extension.")
        ("points-dir", po_value(pov.points_dir), "Points directory.")
        ("points-full-path", po_value(pov.points_full_path), "Full path to the points file. \n (*) Overwrites options: points-name, points-ext, points-dir.")
        ("out-file", po_value(pov.out_file), "Output file.")
        ("resume-after-crash", po_value(pov.resume_after_crash), "Resume after crash.")
        ("check-map-validity", po_value(pov.check_map_validity), "Check map validity.")
        ("algorithm", po_value(pov.algorithm), "Visibility algorithm [simple_polygon, rotational_sweep, triangular_expansion].")
        ("inexact-constructions", po_value(pov.inexact_constructions), "Use inexact constructions.")
        ("max-points", po_value(pov.max_points), "Maximum number of points to test.");
}

char ParseProgramOptions(
    int argc,
    const char *const *argv,
    ProgramOptionVariables &pov
) {
    using namespace trivis_plus::utils;
    po::variables_map vm;
    po::options_description command_line_options;
    po::options_description options_description("Program options");
    AddProgramOptions(options_description, pov);
    try {
        // Parse the command line arguments.
        command_line_options.add(options_description);
        po::store(po::parse_command_line(argc, argv, command_line_options), vm);
        if (vm.count("help")) {
            // If '-h' or '--help' option, print the options and return 'h'.
            command_line_options.print(std::cout, 80);
            return 'h';
        }
        po::notify(vm);
    } catch (const std::exception &e) {
        // If exception, log it and return 'e'.
        InitLogging(severity_level::info);
        LOGF_FTL("Error in parsing arguments: " << e.what() << ".");
        return 'e';
    }
    InitLogging(trivis_plus::utils::severity_level::info);
    // Make file_name, file_extension, file_dir, and file_full_path consistent.
    if (pov.map_full_path.empty()) {
        pov.map_full_path = pov.map_dir + "/" + pov.map_name + (pov.map_extension.empty() ? "" : pov.map_extension);
    } else {
        auto aux = fs::path(pov.map_full_path);
        pov.map_name = aux.replace_extension("").filename().string();
        pov.map_extension = aux.extension().string();
        pov.map_dir = aux.parent_path().string();
    }
    // Make points_name, points_extension, points_dir, and points_full_path consistent.
    if (pov.points_full_path.empty()) {
        pov.points_full_path = pov.points_dir + "/" + pov.points_name + (pov.points_extension.empty() ? "" : pov.points_extension);
    } else {
        auto aux = fs::path(pov.points_full_path);
        pov.points_name = aux.replace_extension("").filename().string();
        pov.points_extension = aux.extension().string();
        pov.points_dir = aux.parent_path().string();
    }
    return '0';
}

template<typename Kernel>
std::pair<std::vector<Segment<Kernel>>, std::vector<std::vector<Segment<Kernel>>>> ToCGAL(const trivis::geom::PolyMap &map) {
    std::vector<Segment<Kernel>> boundary;
    int n_boundary = static_cast<int>(map.border().size());
    boundary.reserve(n_boundary);
    for (int i = 0; i < n_boundary; ++i) {
        boundary.emplace_back(
            Point<Kernel>(map.border()[i].x, map.border()[i].y),
            Point<Kernel>(map.border()[(i + 1) % n_boundary].x, map.border()[(i + 1) % n_boundary].y)
        );
    }
    std::vector<std::vector<Segment<Kernel>>> holes;
    int n_holes = static_cast<int>(map.holes().size());
    holes.reserve(n_holes);
    for (const auto &hole: map.holes()) {
        int n_hole = static_cast<int>(hole.size());
        std::vector<Segment<Kernel>> hole_segments;
        hole_segments.reserve(n_hole);
        for (int i = 0; i < hole.size(); ++i) {
            hole_segments.emplace_back(
                Point<Kernel>(hole[i].x, hole[i].y),
                Point<Kernel>(hole[(i + 1) % n_hole].x, hole[(i + 1) % n_hole].y)
            );
        }
        holes.push_back(std::move(hole_segments));
    }
    return {std::move(boundary), std::move(holes)};
}

template<typename Kernel>
struct LocateQueryResult {
    bool is_inside = false;
    bool is_on_edge = false;
    bool is_on_vertex = false;
    Hh<Kernel> query_halfedge;
};

template<typename Kernel>
LocateQueryResult<Kernel> LocateQuery(
    int id,
    const Point<Kernel> &query,
    const PointLocation<Kernel> &pl,
    Fh<Kernel> interior_face
) {
    LocateQueryResult<Kernel> result;
    auto loc = pl.locate(query);
    if (auto fh = boost::get<Fh<Kernel>>(&loc)) {
        // The query point is in a face.
        if ((*fh) == interior_face) {
            result.is_inside = true;
        }
    } else if (auto eh = boost::get<Hh<Kernel>>(&loc)) {
        result.is_on_edge = true;
        // The query point is on a halfedge.
        if ((*eh)->face() == interior_face) {
            // The halfedge is on the boundary of the interior face.
            result.is_on_edge = true;
            result.query_halfedge = *eh;
        } else if ((*eh)->twin()->face() == interior_face) {
            // The halfedge twin is on the boundary of the interior face.
            result.is_on_edge = true;
            result.query_halfedge = (*eh)->twin();
        }
    } else if (auto vh = boost::get<Vh<Kernel>>(&loc)) {
        result.is_on_vertex = true;
        // The query point is on a vertex.
        if ((*vh)->is_isolated()) {
            if ((*vh)->face() == interior_face) {
                // The vertex is an isolated vertex in the interior face.
                result.is_inside = true;
            }
        } else {
            auto eit = (*vh)->incident_halfedges();
            do {
                if (eit->face() == interior_face) {
                    // The vertex is on the boundary of the interior face.
                    result.is_on_edge = true;
                    result.is_inside = true;
                } else if (eit->twin()->face() == interior_face) {
                    // The vertex is on the boundary of the interior face.
                    eit = eit->twin();
                    result.is_on_edge = true;
                    result.is_inside = true;
                }
                if (result.is_on_edge) {
                    if (eit->target()->point() == query) {
                        result.query_halfedge = eit;
                    } else if (eit->source()->point() == query) {
                        result.query_halfedge = eit->prev();
                    } else {
                        LOGF_ERR("Point " << id << " " << query << " is on a vertex but not on an incident halfedge! This should never happen!");
                    }
                    break;
                }
            } while (++eit != (*vh)->incident_halfedges());
        }
    } else {
        LOGF_ERR("Point " << id << " " << query << " is in an unknown location! This should never happen!");
    }
    return result;
}

template<typename Kernel>
struct VisibilityQueryResult {
    Fh<Kernel> fh;
    Arrangement<Kernel> visi_polygon;
};

template<typename Kernel, typename Visibility>
int MainBody(const ProgramOptionVariables &pov) {

    trivis::utils::SimpleClock clock;
    std::stringstream info;

    clock.Restart();
    auto map_opt = trivis_plus::data_loading::LoadPolyMap(pov.map_full_path, std::nullopt, &info);
    if (!map_opt) {
        LOGF_FTL("Error while loading map: '" << info.str() << "'.");
        return EXIT_FAILURE;
    }
    auto map = std::move(map_opt.value());
    map_opt = std::nullopt; // cannot use map_opt anymore ! (it was moved)
    std::ifstream ifs(pov.points_full_path);
    if (!ifs.is_open()) {
        LOGF_FTL("Cannot open file " << pov.points_full_path << " for reading.");
        return EXIT_FAILURE;
    }

    int n_points;
    ifs >> n_points;
    n_points = std::min(n_points, pov.max_points);

    std::ofstream ofs;
    int start_id = 0;
    if (pov.resume_after_crash) {
        {   // Read the last point id from the output file.
            std::ifstream ifs_out(pov.out_file);
            if (!ifs_out.is_open()) {
                LOGF_FTL("Cannot open file " << pov.out_file << " for reading.");
                return EXIT_FAILURE;
            }
            std::string line;
            while (ifs_out >> std::ws && std::getline(ifs_out, line));
            std::istringstream{line} >> start_id;
            ++start_id; // The next point id.
        }
        ofs = std::ofstream(pov.out_file, std::ios_base::app);
        if (!ofs.is_open()) {
            LOGF_FTL("Cannot open file " << pov.out_file << " for appending.");
            return EXIT_FAILURE;
        }
    } else {
        ofs = std::ofstream(pov.out_file);
        if (!ofs.is_open()) {
            LOGF_FTL("Cannot open file " << pov.out_file << " for writing.");
            return EXIT_FAILURE;
        }
    }

    auto [boundary, holes] = ToCGAL<Kernel>(map);

    clock.Restart();
    Arrangement<Kernel> env;
    CGAL::insert_non_intersecting_curves(env, boundary.begin(), boundary.end());
    for (const auto &hole: holes) {
        CGAL::insert_non_intersecting_curves(env, hole.begin(), hole.end());
    }
    Fh<Kernel> interior_face = *(++env.face_handles().begin());
    double time_init = clock.TimeInSeconds();

    if (!pov.resume_after_crash) {
        ofs << "time_init" << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_init << "\n";
    }
    if (pov.check_map_validity) {
        if (!env.is_valid()) {
            LOGF_WRN("Arrangement is not valid!");
        }
    }

    clock.Restart();
    PointLocation<Kernel> pl(env);
    Visibility vis(env);
    double time_preprocess = clock.TimeInSeconds();

    if (!pov.resume_after_crash) {
        ofs << "time_preprocess" << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_preprocess << "\n";
        ofs << n_points << "\n";
        ofs.flush();
    }

    int id;
    double x, y;
    for (int i = 0; i < n_points; ++i) {
        ifs >> id >> x >> y;
        if (id < start_id) {
            continue;
        }
        clock.Restart();
        auto query = Point<Kernel>(x, y);
        auto locate_result = LocateQuery(id, query, pl, interior_face);
        double time_is_inside = clock.TimeInSeconds();
        clock.Restart();
        std::optional<VisibilityQueryResult<Kernel>> visi_result = std::nullopt;
        if (locate_result.is_inside) {
            visi_result = VisibilityQueryResult<Kernel>();
            if (locate_result.is_on_edge) {
                if (locate_result.is_on_vertex && locate_result.query_halfedge->target()->point() != query) {
                    LOGF_ERR("Point " << id << " is on a vertex but the target of the halfedge is not the query point! This should never happen!");
                }
                if (locate_result.query_halfedge->face() != interior_face) {
                    LOGF_ERR("Point " << id << " is on an edge but the face of the halfedge is not the interior face! This should never happen!");
                }
                visi_result->fh = vis.compute_visibility(query, locate_result.query_halfedge, visi_result->visi_polygon);
            } else {
                visi_result->fh = vis.compute_visibility(query, interior_face, visi_result->visi_polygon);
            }
        }
        double time_vis_poly = clock.TimeInSeconds();
        ofs << id << " times " << 2
            << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_is_inside
            << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_vis_poly
            << " polygon ";
        if (visi_result) {
            std::vector<Point<Kernel>> polygon;
            auto eit = visi_result->fh->outer_ccb();
            do {
                polygon.emplace_back(eit->target()->point());
            } while (++eit != visi_result->fh->outer_ccb());
            ofs << polygon.size();
            for (const auto &p: polygon) {
                ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << p.x();
                ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << p.y();
            }
        } else {
            ofs << 0;
        }
        ofs << "\n";
        ofs.flush();
    }
    ofs << "DONE\n";
    ofs.flush();
    ofs.close();

    return EXIT_SUCCESS;
}

int main(
    int argc,
    const char *const *argv
) {
    ProgramOptionVariables pov;
    char c = ParseProgramOptions(argc, argv, pov);
    if (c == 'h') return EXIT_SUCCESS;
    if (c == 'e') return EXIT_FAILURE;
    if (pov.inexact_constructions) {
        using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
        if (pov.algorithm == "simple_polygon") return MainBody<Kernel, VisibilitySimplePolygon<Kernel>>(pov);
        if (pov.algorithm == "rotational_sweep") return MainBody<Kernel, VisibilityRotationalSweep<Kernel>>(pov);
        if (pov.algorithm == "triangular_expansion") return MainBody<Kernel, VisibilityTriangularExpansion<Kernel>>(pov);
    } else {
        using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
        if (pov.algorithm == "simple_polygon") return MainBody<Kernel, VisibilitySimplePolygon<Kernel>>(pov);
        if (pov.algorithm == "rotational_sweep") return MainBody<Kernel, VisibilityRotationalSweep<Kernel>>(pov);
        if (pov.algorithm == "triangular_expansion") return MainBody<Kernel, VisibilityTriangularExpansion<Kernel>>(pov);
    }
    LOGF_FTL("Unknown visibility algorithm: " << pov.algorithm);
    return EXIT_FAILURE;
}