/**
 * File:   generate_meshes.cc
 *
 * Date:   03.10.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include <iomanip>
#include <fstream>

// === BOOST INCLUDES ===

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// === TRIVIS INCLUDES ===

#include "trivis/trivis.h"
#include "trivis/utils/simple_clock.h"
#include "trivis/utils/clipper_utils.h"
#include "trivis/utils/random_points.h"
#include "trivis/geom/robust_geometry.h"
#include "trivis/geom/generic_geom_types.h"
#include "trivis/geom/generic_geom_utils.h"
#include "trivis/geom/intersections.h"

#include "trivis_plus/data_loading/load_map.h"
#include "trivis_plus/data_loading/load_mesh.h"
#include "trivis_plus/drawing/drawing.h"
#include "trivis_plus/drawing/random_colors.h"
#include "trivis_plus/utils/log.h"
#include "trivis_plus/utils/ptree_json.h"

#include "mesh_optim/mesh_weights.h"
#include "mesh_optim/min_weight_triangulation.h"
#include "mesh_optim/triangular_meshes.h"

// === THIS PROJECT INCLUDES ===

#ifndef DEFAULT_MAP_DIR
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_OUT_DIR
#define DEFAULT_OUT_DIR "."
#endif

#ifndef DEFAULT_WEIGHTS_DIR
#define DEFAULT_WEIGHTS_DIR "."
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace dr = trivis_plus::drawing;

/**
 * All program option variables and their default values should be defined here.
 * For each variable, there should be an option added in AddProgramOptions.
 */
struct ProgramOptionVariables {
    bool verbosity_silent = false;
    bool verbosity_quiet = false;
    bool verbosity_verbose = false;
    std::string mesh_file;
    std::string mesh_info_file;
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string output_dir = DEFAULT_OUT_DIR;
    double bucket_size = 0.75;
    double mean_bucket_triangle_count_max = -1.0;
    int n_queries = 1000000;
    unsigned random_seed_queries = 42;
    bool evaluate_tea = false;
    bool evaluate_tea_early_exit = false;
    bool evaluate_edges = false;
    std::vector<double> vis_radii_test = {-1.0};
    bool match_radii = false;
    bool match_radii_ignore_inf = false;
    double expansion_estimate_sample_dist = -1.0;
    bool save_avg_expansions = false;
    bool save_all_expansions = false;
};

template<typename T>
std::string Vec2Str(const std::vector<T> &vec) {
    std::stringstream ret;
    for (int i = 0; i < vec.size(); ++i) {
        ret << vec[i] << ((i == vec.size() - 1) ? "" : ", ");
    }
    return ret.str();
}

void AddProgramOptions(
    po::options_description &options_description,
    ProgramOptionVariables &pov
) {
    options_description.add_options()
        ("help,h", "Produce this help message. \n (*) Overwrites options: all.")
        ("verbose",
         po::bool_switch(&pov.verbosity_verbose)->default_value(pov.verbosity_verbose),
         "Log messages: fatal, error, warning, info, debug.")
        ("quiet",
         po::bool_switch(&pov.verbosity_quiet)->default_value(pov.verbosity_quiet),
         "Log messages: fatal, error, warning. \n (*) Overwrites options: verbose.")
        ("silent",
         po::bool_switch(&pov.verbosity_silent)->default_value(pov.verbosity_silent),
         "Log messages: fatal, error. \n (*) Overwrites options: verbose, quiet.")
        ("mesh-file",
         po::value(&pov.mesh_file)->required(),
         "[REQUIRED] TXT file that contains the mesh.")
        ("mesh-info-file",
         po::value(&pov.mesh_info_file)->required(),
         "[REQUIRED] JSON file that contains info about the mesh.")
        ("map-ext",
         po::value(&pov.map_extension)->default_value(pov.map_extension),
         "Map file extension.")
        ("map-dir",
         po::value(&pov.map_dir)->default_value(pov.map_dir),
         "Map file directory.")
        ("out-dir",
         po::value(&pov.output_dir)->default_value(pov.output_dir),
         "Output directory.")
        ("bucket-size",
         po::value(&pov.bucket_size)->default_value(pov.bucket_size),
         "Size of the buckets used for locating a point in the triangular mesh.")
        ("mean-bucket-triangle-count-max",
         po::value(&pov.mean_bucket_triangle_count_max),
         "The maximum allowed mean value of how many triangle should be in each bucket. \n (*) Overwrites options: bucket-size.")
        ("n-queries",
         po::value(&pov.n_queries)->default_value(pov.n_queries),
         "How many random query points should be generated.")
        ("random-seed-queries",
         po::value(&pov.random_seed_queries)->default_value(pov.random_seed_queries),
         "Seed for the random generator of query points.")
        ("evaluate-tea",
         po::bool_switch(&pov.evaluate_tea)->default_value(pov.evaluate_tea),
         "This program will evaluate the TEA.")
        ("evaluate-edges",
         po::bool_switch(&pov.evaluate_edges)->default_value(pov.evaluate_edges),
         "This program will evaluate edges of the mesh. \n (*) Overwrites options: ")
        ("evaluate-tea-early-exit",
         po::bool_switch(&pov.evaluate_tea_early_exit)->default_value(pov.evaluate_tea_early_exit),
         "This program will evaluate the TEA with early exit (TEA+).")
        ("vis-radii-test",
         po::value(&pov.vis_radii_test)->multitoken()->default_value(pov.vis_radii_test, "[" + Vec2Str(pov.vis_radii_test) + "]"),
         "Visibility radii used for testing the triangulations (-1 ~ infinite).")
        ("match-radii",
         po::bool_switch(&pov.match_radii)->default_value(pov.match_radii),
         "Methods starting with d- will only be tested on computing vis. regions with the same limited vis. radius d.")
        ("match-radii-ignore-inf",
         po::bool_switch(&pov.match_radii_ignore_inf)->default_value(pov.match_radii_ignore_inf),
         "Methods starting with d- will only be tested on computing vis. regions with the same limited vis. radius d, and d=inf.")
        ("compute-expansion-estimate",
         po::value(&pov.expansion_estimate_sample_dist)->default_value(pov.expansion_estimate_sample_dist),
         "If this value is set to >0, this program computes the precise estimate of the mean number of expanded edges using the value as sampling distance and then terminates.")
        ("save-avg-expansions",
         po::bool_switch(&pov.save_avg_expansions)->default_value(pov.save_avg_expansions),
         "If used, average data from each query are saved to a special file (avg_expansions_(...).txt).")
        ("save-all-expansions",
         po::bool_switch(&pov.save_all_expansions)->default_value(pov.save_all_expansions),
         "If used, all data from each query are saved to a special file (all_expansions_(...).txt).");
}

trivis_plus::utils::severity_level GetSeverity(const ProgramOptionVariables &pov) {
    using namespace trivis_plus::utils;
    if (pov.verbosity_silent) {
        return severity_level::error;
    } else if (pov.verbosity_quiet) {
        return severity_level::warning;
    } else if (pov.verbosity_verbose) {
        return severity_level::debug;
    } else {
        return severity_level::info;
    }
}

/**
 *
 * Parses arguments and initializes logging.
 *
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @param pov
 * @return Character 'e' if an exception occurred, 'h' if --help option, '0' else.
 */
char ParseProgramOptions(
    int argc,
    const char *const *argv,
    ProgramOptionVariables &pov
) {
    using namespace trivis_plus::utils;
    if (argc < 2) {
        InitLogging(severity_level::info);
        LOGF_FTL("No arguments provided. Use --help to see available options.");
        return 'e';
    }
    po::variables_map vm;
    po::options_description command_line_options;
    po::options_description options_description("General options");
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
    InitLogging(GetSeverity(pov));
    return '0';
}

int EvaluateEdges(
    const ProgramOptionVariables &pov,
    const trivis::Trivis &vis
) {

    std::vector<double> dist_ratios;
    double longest_possible_non_obstacle_edge_length = 0.0;
    for (int node_id = 0; node_id < vis.mesh().vertices.size(); ++node_id) {
        const auto &node = vis.mesh().vertices[node_id];
        std::vector<int> obstacle_edges;
        for (int e_id: node.edges) {
            if (vis.mesh().edges[e_id].is_boundary()) {
                obstacle_edges.push_back(e_id);
            }
        }
        std::vector<int> visible_vertices = vis.VisibleVertices(node_id);
        for (int visible_node_id: visible_vertices) {
            bool is_boundary = false;
            for (int e_id: obstacle_edges) {
                const auto &e = vis.mesh().edges[e_id];
                if (e.vertices[e.vertices[0] == node_id ? 1 : 0] == visible_node_id) {
                    is_boundary = true;
                    break;
                }
            }
            if (!is_boundary) {
                double dist = node.point.DistanceTo(vis.mesh().point(visible_node_id));
                if (dist > longest_possible_non_obstacle_edge_length) {
                    longest_possible_non_obstacle_edge_length = dist;
                }
            }
        }
    }
    for (int e_id = 0; e_id < vis.mesh().edges.size(); ++e_id) {
        const auto &e = vis.mesh().edges[e_id];
        if (e.is_boundary()) {
            continue;
        }
        double e_length = vis.mesh().point(e.vertices[0]).DistanceTo(vis.mesh().point(e.vertices[1]));
        double dist_ratio = e_length / longest_possible_non_obstacle_edge_length;
        dist_ratios.push_back(dist_ratio);
    }

    double min = *std::min_element(dist_ratios.begin(), dist_ratios.end());
    double max = *std::max_element(dist_ratios.begin(), dist_ratios.end());
    double mean = std::reduce(dist_ratios.begin(), dist_ratios.end()) / static_cast<double>(dist_ratios.size());
    LOGF_INF("Min: " << min << ", mean: " << mean << ", max: " << max << ".");

    return EXIT_SUCCESS;
}

double ComputeExpansionEstimate(
    const ProgramOptionVariables &pov,
    const trivis::Trivis &vis
) {
    double ret = 0.0;
    double total_area = Clipper2Lib::Area(Clipper2Lib::PathsDToPaths64(trivis::utils::Geom2Clipper(vis.triangles())));
    for (int edge_id = 0; edge_id < vis.mesh().edges.size(); ++edge_id) {
        const auto &e = vis.mesh().edges[edge_id];
        if (e.is_boundary()) {
            continue;
        }
        const auto &p0 = vis.mesh().point(e.vertices[0]);
        const auto &p1 = vis.mesh().point(e.vertices[1]);
        auto samples = trivis::geom::SampleSegment2D(p0, p1, pov.expansion_estimate_sample_dist, false, false);
        if (samples.empty()) {
            samples.push_back((p0 + p1) / 2.0);
        }
        Clipper2Lib::Paths64 solution;
        Clipper2Lib::Paths64 vis_polygons64;
        Clipper2Lib::Clipper64 clipper64;
        for (int i = 0; i < samples.size(); ++i) {
            const auto &s = samples[i];
            auto s_locate = vis.LocatePoint(s);
            if (!s_locate) {
                continue;
            }
            auto abs_reg_opt = vis.VisibilityRegion(s, s_locate.value());
            vis_polygons64.push_back(Clipper2Lib::PathDToPath64(trivis::utils::Geom2Clipper(vis.ToRadialVisibilityRegion(abs_reg_opt).ToPolygon())));
            if (vis_polygons64.size() == 2 || i + 1 == static_cast<int>(samples.size())) {
                clipper64.AddClip(solution);
                clipper64.AddSubject(vis_polygons64);
                clipper64.Execute(Clipper2Lib::ClipType::Union, Clipper2Lib::FillRule::NonZero, solution);
                clipper64.Clear();
                vis_polygons64.clear();
            }
        }
        ret += Clipper2Lib::Area(solution) / total_area;
    }
    return ret;
}

/**
 *
 *  ##########################################
 *  ## THIS IS THE MAIN BODY OF THE PROGRAM ##
 *  ##########################################
 *
 * @param pov ~ program option variables
 * @return exit code
 */
int MainBody(const ProgramOptionVariables &pov) {

    fs::create_directories(pov.output_dir);

    trivis::utils::SimpleClock clock;

    LOGF_INF("Loading the mesh.");
    // Load the mesh.
    trivis::mesh::TriMesh mesh;
    {
        std::stringstream info_ss;
        auto mesh_opt = trivis_plus::data_loading::LoadTriMesh(pov.mesh_file, &info_ss);
        if (!mesh_opt) {
            LOGF_FTL("Could not load the mesh from file " << pov.mesh_file << ". Message: " << info_ss.str());
            return EXIT_FAILURE;
        }
        mesh = std::move(*mesh_opt);
    }

    LOGF_INF("Loading the mesh info.");
    // Load some information about the mesh.
    std::string map_name;
    double map_scale;
    std::string method_name;
    double method_vis_radius;
    trivis_plus::utils::PTreeJSON mesh_info;
    try {
        mesh_info.Read(pov.mesh_info_file);
        mesh_info.Get("map_name", map_name);
        mesh_info.Get("map_scale", map_scale);
        mesh_info.Get("method_name", method_name);
        mesh_info.Get("method_vis_radius", method_vis_radius);
    } catch (const std::exception &e) {
        LOGF_FTL("Could not read file" << pov.mesh_info_file << ".");
        return EXIT_FAILURE;
    }

    LOGF_INF("Initializing Trivis.");
    // Create and initialize the Trivis object.
    trivis::Trivis vis;
    {   // Load map from file and move it to Trivis (without copying).
        std::stringstream info;
        auto map = trivis_plus::data_loading::LoadPolyMap(pov.map_dir + "/" + map_name + pov.map_extension, map_scale, &info);
        if (!map.has_value()) {
            LOGF_FTL("Error while loading map. " << info.str());
            return EXIT_FAILURE;
        }
        // Warning: order of the following matters!
        map->ShiftToOrigin(); // Subtracts min X and Y coordinates from all points.
        map->RemoveDuplicatePoints(); // Removes all consecutive identical points.
        map->SimplifyWeaklySelfIntersectingPolygons(); // Splits all weakly simple polygons to multiple (touching) strongly simple polygons.
        map->RemoveCollinearPoints(); // Removes all consecutive collinear points.
        vis.SetMap(std::move(map.value()));
        // cannot use map anymore ! (it was moved)
    }
    vis.ConstructMeshCDT();
    auto cdt_triangles = vis.triangles();

    // Set the mesh to Trivis.
    LOGF_INF("Setting the mesh to Trivis.");
    clock.Restart();
    vis.SetMesh(mesh);
    if (pov.mean_bucket_triangle_count_max > 0.0) {
        vis.FillPointLocationBuckets(std::nullopt, pov.mean_bucket_triangle_count_max);
    } else {
        vis.FillPointLocationBuckets(pov.bucket_size);
    }
    vis.OptimizePointLocationBucketTriangles();
    double time_vis_init = clock.TimeInSeconds();

    if (pov.expansion_estimate_sample_dist > 0.0) {
        LOGF_INF("Computing expansions estimate.");
        double expansion_estimate = ComputeExpansionEstimate(pov, vis);
        auto out_file = pov.output_dir + "/exp_estimates.txt";
        std::ofstream ofs;
        ofs.open(out_file, std::ios_base::app);
        ofs << std::setprecision(std::numeric_limits<double>::max_digits10) << pov.expansion_estimate_sample_dist;
        ofs << " " << std::setprecision(std::numeric_limits<double>::max_digits10) << expansion_estimate;
        ofs << "\n";
        ofs.close();
        LOGF_INF("Expansions estimate: " << expansion_estimate);
        return EXIT_SUCCESS;
    }

    // Evaluate the mesh.
    LOGF_INF("Evaluating the mesh.");
    if (pov.evaluate_edges) {
        return EvaluateEdges(pov, vis);
    }
    std::vector<bool> early_exits;
    if (pov.evaluate_tea) early_exits.push_back(false);
    if (pov.evaluate_tea_early_exit) early_exits.push_back(true);
    if (early_exits.empty()) {
        LOGF_WRN("Nothing to evaluate.");
        return EXIT_SUCCESS;
    }

    LOGF_INF("Precomputing " << pov.n_queries << " random query points.");
    // Precompute the random points.
    trivis::geom::FPoints random_points;
    {
        std::vector<double> accum_areas;
        auto rng_points = std::mt19937(pov.random_seed_queries);
        random_points.reserve(pov.n_queries);
        for (int i = 0; i < pov.n_queries; ++i) {
            random_points.push_back(trivis::utils::UniformRandomPointInRandomTriangle(cdt_triangles, accum_areas, rng_points));
        }
    }

    for (bool early_exit: early_exits) {
        if (early_exit) {
            LOGF_INF("Evaluating d-TEA (with early exit).");
        } else {
            LOGF_INF("Evaluating TEA.");
        }
        for (double vis_radius: pov.vis_radii_test) {

            if ((pov.match_radii || pov.match_radii_ignore_inf) && method_vis_radius >= 0.0 &&
                method_vis_radius != vis_radius) {
                if (!pov.match_radii_ignore_inf || vis_radius > 0.0) {
                    continue;
                }
            }

            LOGF_INF("Evaluating d = " << (vis_radius > 0.0 ? std::to_string(vis_radius) : "inf") << ".");

            trivis_plus::utils::PTreeJSON records;
            records.Put("vis_radius", vis_radius);
            records.Put("n_queries", pov.n_queries);
            records.Put("random_seed_queries", pov.random_seed_queries);
            records.Put("method_early_exit", early_exit);
            records.Put("time_vis_init", time_vis_init);
            records.Put("bucket_size", pov.bucket_size);
            records.Put("mean_bucket_triangle_count_max", pov.mean_bucket_triangle_count_max);

            int n_r = static_cast<int>(random_points.size());

            int num_expansions_total = 0;
            double time_vis_reg_bucketing = 0.0;
            double time_vis_reg_tea = 0.0;
            double time_vis_reg_intersections = 0.0;
            double time_vis_reg_radial = 0.0;

            std::optional<std::ofstream> avg_expansions_file_opt;
            if (pov.save_avg_expansions) {
                std::string avg_expansions_file_str =
                    pov.output_dir + "/avg_expansions_" + (early_exit ? "d" : "") + "tea_d-" + (vis_radius > 0.0 ? std::to_string(vis_radius) : "inf") + ".txt";
                avg_expansions_file_opt = std::ofstream{};
                avg_expansions_file_opt->open(avg_expansions_file_str);
            }

            std::optional<std::ofstream> all_expansions_file_opt;
            if (pov.save_all_expansions) {
                std::string all_expansions_file_str =
                    pov.output_dir + "/all_expansions_" + (early_exit ? "d" : "") + "tea_d-" + (vis_radius > 0.0 ? std::to_string(vis_radius) : "inf") + ".txt";
                all_expansions_file_opt = std::ofstream{};
                all_expansions_file_opt->open(all_expansions_file_str);
            }

            trivis::utils::SimpleClock clock_total;
            int cnt = 0;
            for (const auto &p: random_points) {
                ++cnt;
                clock.Restart();
                auto p_locate = vis.LocatePoint(p);
                double time_vis_reg_bucketing_curr = clock.TimeInSeconds();
                time_vis_reg_bucketing += time_vis_reg_bucketing_curr;
                // Compute an abstract visibility region (a structure without computed intersections).
                trivis::Trivis::ExpansionStats stats;
                clock.Restart();
                if (!p_locate) {
                    LOGF_WRN("Could not compute visibility region for query point " << p << ".");
                    continue;
                }
                trivis::AbstractVisibilityRegion reg_abstract;
                if (early_exit && vis_radius > 0.0) {
                    reg_abstract = vis.VisibilityRegion(p, p_locate.value(), vis_radius, &stats);
                } else {
                    reg_abstract = vis.VisibilityRegion(p, p_locate.value(), std::nullopt, &stats);
                }
                double time_vis_reg_tea_curr = clock.TimeInSeconds();
                time_vis_reg_tea += time_vis_reg_tea_curr;
                num_expansions_total += stats.num_expansions;
                // Convert the abstract to a more concrete representation (compute the intersections).
                clock.Restart();
                auto reg = vis.ToRadialVisibilityRegion(reg_abstract);
                double time_vis_reg_intersections_curr = clock.TimeInSeconds();
                time_vis_reg_intersections += time_vis_reg_intersections_curr;
                // Convert the concrete to proper radial region (intersect it with a circle of radius vis_radius).
                double time_vis_reg_radial_curr = 0.0;
                if (vis_radius > 0.0) {
                    clock.Restart();
                    reg.IntersectWithCircleCenteredAtSeed();
                    time_vis_reg_radial_curr = clock.TimeInSeconds();
                    time_vis_reg_radial += time_vis_reg_radial_curr;
                }
                if (avg_expansions_file_opt) {
                    *avg_expansions_file_opt << std::setprecision(std::numeric_limits<double>::max_digits10) << static_cast<double>(num_expansions_total) / cnt << "\n";
                }
                if (all_expansions_file_opt) {
                    *all_expansions_file_opt << stats.num_expansions;
                    *all_expansions_file_opt << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_vis_reg_bucketing_curr;
                    *all_expansions_file_opt << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_vis_reg_tea_curr;
                    *all_expansions_file_opt << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_vis_reg_intersections_curr;
                    *all_expansions_file_opt << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_vis_reg_radial_curr;
                    *all_expansions_file_opt << "\n";
                }
            }
            if (avg_expansions_file_opt) {
                avg_expansions_file_opt->close();
            }
            if (all_expansions_file_opt) {
                all_expansions_file_opt->close();
            }
            double time_vis_reg_total = clock_total.TimeInSeconds();
            double avg_expansions = static_cast<double>(num_expansions_total) / n_r;
            LOGF_INF("Avg. expansions: " << avg_expansions);
            records.Put("avg_expansions", avg_expansions);
            records.Put("avg_time_vis_reg_bucketing", time_vis_reg_bucketing / n_r);
            records.Put("avg_time_vis_reg_tea", time_vis_reg_tea / n_r);
            records.Put("avg_time_vis_reg_intersections", time_vis_reg_intersections / n_r);
            records.Put("avg_time_vis_reg_radial", time_vis_reg_radial / n_r);
            records.Put("avg_time_vis_reg_total", time_vis_reg_total / n_r);

            std::string record_file = pov.output_dir + "/mesh_eval_" + (early_exit ? "d" : "") + "tea_d-" + (vis_radius > 0.0 ? std::to_string(vis_radius) : "inf") + ".json";
            try {
                records.Write(record_file);
                LOGF_INF("Records written to " << record_file << ".");
            } catch (const std::exception &e) {
                LOGF_ERR("Could not write to file " << record_file << ".");
            }

        }
    }

    LOGF_INF("DONE");

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
    return MainBody(pov);
}