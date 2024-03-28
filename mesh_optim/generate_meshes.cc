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
    std::string map_name = "undefined";
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string output_dir = DEFAULT_OUT_DIR;
    std::string weights_dir = DEFAULT_WEIGHTS_DIR;
    bool print_map_info = false;
    bool draw = false;
    bool save_weights = false;
    bool load_weights = false;
    bool precompute_weights = false;
    bool dont_save_mesh = false;
    bool save_only_last_iteration = false;
    double map_scale = 0.01;
    double bucket_size = 0.75;
    double mean_bucket_triangle_count_max = -1.0;
    std::vector<std::string> methods = {"CDT", "MinLT", "MaxLT", "MinVT", "MaxVT", "d-MinVT", "d-MaxVT"};
    std::vector<double> vis_radii_mesh;
    double max_sample_dist = 0.1;
    int mwt_iter_limit = 10;
    double mwt_improvement_ratio_limit = -1.0;
    double mwt_time_limit = -1.0;
    int mwt_simple_polygon_max_size = -1;
    unsigned random_seed_method = std::random_device{}();
    double eps_dist_diff_collinear = 1e-7;
    double weights_long_edge_penalty = -1.0;
    double weights_long_edge_threshold = -1.0;
    int weights_shortest_edges_ver_max = -1;
};

template<typename T>
std::string Vec2Str(const std::vector<T> &vec) {
    std::stringstream ret;
    for (int i = 0; i < vec.size(); ++i) {
        ret << vec[i] << ((i == vec.size() - 1) ? "" : " ");
    }
    return ret.str();
}

struct MeshWithInfo {
    std::string name;
    trivis::mesh::TriMesh mesh;
    double time_weights = 0.0;
    double time_construction = 0.0;
    std::optional<long> random_seed = std::nullopt;
    std::optional<double> weight = std::nullopt;
    std::optional<double> vis_radius = std::nullopt;
    std::optional<int> iteration = std::nullopt;
};

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
        ("map-name",
         po::value(&pov.map_name)->default_value(pov.map_name),
         "Map name.")
        ("map-ext",
         po::value(&pov.map_extension)->default_value(pov.map_extension),
         "Map file extension.")
        ("map-dir",
         po::value(&pov.map_dir)->default_value(pov.map_dir),
         "Map file directory.")
        ("map",
         po::value(&pov.map_full_path)->default_value(pov.map_full_path),
         "Full path to the map file. \n (*) Overwrites options: map-name, map-ext, map-dir.")
        ("out-dir",
         po::value(&pov.output_dir)->default_value(pov.output_dir),
         "Output directory.")
        ("weights-dir",
         po::value(&pov.weights_dir)->default_value(pov.weights_dir),
         "Directory for storing precomputed weights.")
        ("print-map-info",
         po::bool_switch(&pov.print_map_info)->default_value(pov.print_map_info),
         "Prints map info, then the program is terminated.")
        ("draw",
         po::bool_switch(&pov.draw)->default_value(pov.draw),
         "Enables drawing of the map and mesh.")
        ("save-weights",
         po::bool_switch(&pov.save_weights)->default_value(pov.save_weights),
         "The weights will be saved to weights_dir directory.")
        ("load-weights",
         po::bool_switch(&pov.load_weights)->default_value(pov.load_weights),
         "Instead of computing the weights, the weights will be loaded from weights_dir directory.")
        ("precompute-weights",
         po::bool_switch(&pov.precompute_weights)->default_value(pov.precompute_weights),
         "The weights will always be computed and saved to weights_dir directory, then the program is terminated. \n (*) Overwrites options: save-weights, load-weights.")
        ("dont-save-mesh",
         po::bool_switch(&pov.dont_save_mesh)->default_value(pov.dont_save_mesh),
         "The mesh won't be saved to mesh.txt file. File mesh_info.json will still be generated.")
        ("save-only-last-iteration",
         po::bool_switch(&pov.save_only_last_iteration)->default_value(pov.save_only_last_iteration),
         "Only the last iteration is saved.")
        ("map-scale",
         po::value(&pov.map_scale),
         "Map scale (optional).")
        ("bucket-size",
         po::value(&pov.bucket_size)->default_value(pov.bucket_size),
         "Size of the buckets used for locating a point in the triangular mesh.")
        ("mean-bucket-triangle-count-max",
         po::value(&pov.mean_bucket_triangle_count_max),
         "The maximum allowed mean value of how many triangle should be in each bucketmak. \n (*) Overwrites options: bucket-size.")
        ("methods",
         po::value(&pov.methods)->multitoken()->default_value(pov.methods, Vec2Str(pov.methods)),
         "Methods used in the tests.")
        ("vis-radii-mesh",
         po::value(&pov.vis_radii_mesh)->multitoken()->default_value(pov.vis_radii_mesh, Vec2Str(pov.vis_radii_mesh)),
         "Visibility radii used for constructing the triangulations (-1 ~ infinite).")
        ("max-sample-dist",
         po::value(&pov.max_sample_dist)->default_value(pov.max_sample_dist),
         "Max. sample dist for computing visibility regions from segments.")
        ("mwt-iter-limit",
         po::value(&pov.mwt_iter_limit)->default_value(pov.mwt_iter_limit),
         "MWT construction no. of iterations limit.")
        ("mwt-improvement-ratio-limit",
         po::value(&pov.mwt_improvement_ratio_limit)->default_value(pov.mwt_improvement_ratio_limit),
         "MWT construction improvement ratio limit.")
        ("mwt-time-limit",
         po::value(&pov.mwt_time_limit)->default_value(pov.mwt_time_limit),
         "MWT construction time limit.")
        ("mwt-simple-polygon-max-size",
         po::value(&pov.mwt_simple_polygon_max_size)->default_value(pov.mwt_simple_polygon_max_size),
         "Maximum size of the random simple polygons generated when computing MWT.")
        ("random-seed-method",
         po::value(&pov.random_seed_method)->default_value(pov.random_seed_method),
         "Seed for the random generator inside the method (generated by std::random_device{}() by default).")
        ("eps-dist-diff-collinear",
         po::value(&pov.eps_dist_diff_collinear)->default_value(pov.eps_dist_diff_collinear),
         "Epsilon for determining close-to-collinear points based on their pair-wise distance.")
        ("weights-long-edge-penalty",
         po::value(&pov.weights_long_edge_penalty)->default_value(pov.weights_long_edge_penalty),
         "Positive value (-1 ~ unset). If this value is set, too long edges will by penalized by assigning this value to their weight instead of computing it standardly.")
        ("weights-long-edge-threshold",
         po::value(&pov.weights_long_edge_threshold)->default_value(pov.weights_long_edge_threshold),
         "Value in range [0, 1]. 100 times this value % of edges are considered ok while the rest is too long.")
        ("weights-shortest-edges-ver-max",
         po::value(&pov.weights_shortest_edges_ver_max)->default_value(pov.weights_shortest_edges_ver_max),
         "An experimental parameter affecting weights computation. Do not use.");
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
    // Make file_name, file_extension, file_dir, and file_full_path consistent.
    if (pov.map_full_path.empty()) {
        pov.map_full_path = pov.map_dir + "/" + pov.map_name + (pov.map_extension.empty() ? "" : pov.map_extension);
    } else {
        auto aux = fs::path(pov.map_full_path);
        pov.map_name = aux.replace_extension("").filename().string();
        pov.map_extension = aux.extension().string();
        pov.map_dir = aux.parent_path().string();
    }
    return '0';
}

inline bool Contains(
    const std::set<std::string> &m_set,
    const std::string &entry
) {
    return m_set.find(entry) != m_set.end();
}

void SaveMeshWithInfo(
    const ProgramOptionVariables &pov,
    const MeshWithInfo &mwi,
    const std::string &output_dir,
    dr::MapDrawer *drawer = nullptr
) {
    fs::create_directories(output_dir);

    if (drawer) {
        drawer->OpenPDF(output_dir + "/mesh.pdf");
        drawer->DrawMap(dr::kColorLightSkyBlue, dr::kColorWhite);
        drawer->DrawPolygons(Mesh2Polygons(mwi.mesh), 0.5, dr::kColorBlack);
        drawer->Close();
    }

    if (!pov.dont_save_mesh) {
        std::stringstream info_ss;
        bool ok = trivis_plus::data_loading::SaveTriMesh(mwi.mesh, output_dir + "/mesh.txt", &info_ss);
        if (!ok) {
            LOGF_WRN("Could not save mesh. Info: " << info_ss.str());
        }
    }

    trivis_plus::utils::PTreeJSON info;
    info.Put("map_name", pov.map_name);
    info.Put("map_scale", pov.map_scale);
    info.Put("max_sample_dist", pov.max_sample_dist);
    info.Put("method_name", mwi.name);
    info.Put("method_vis_radius", mwi.vis_radius ? static_cast<double>(*mwi.vis_radius) : -1.0);
    info.Put("method_random_seed", mwi.random_seed ? static_cast<int>(*mwi.random_seed) : -1);
    info.Put("weights_bucket_size", pov.bucket_size);
    info.Put("weights_mean_bucket_triangle_count_max", pov.mean_bucket_triangle_count_max);
    info.Put("weights_long_edge_penalty", pov.weights_long_edge_penalty);
    info.Put("weights_long_edge_threshold", pov.weights_long_edge_threshold);
    info.Put("weights_shortest_edges_ver_max", pov.weights_shortest_edges_ver_max);
    info.Put("mwt_iteration", mwi.iteration ? *mwi.iteration : 1);
    info.Put("mwt_iter_limit", pov.mwt_iter_limit);
    info.Put("mwt_improvement_ratio_limit", pov.mwt_improvement_ratio_limit);
    info.Put("mwt_time_limit", pov.mwt_time_limit);
    info.Put("mwt_simple_polygon_max_size", pov.mwt_simple_polygon_max_size);
    info.Put("mwt_result_weight", mwi.weight ? *mwi.weight : std::numeric_limits<double>::quiet_NaN());
    info.Put("eps_dist_diff_collinear", pov.eps_dist_diff_collinear);
    info.Put("time_weights", mwi.time_weights);
    info.Put("time_construction", mwi.time_construction);

    std::string record_file = output_dir + "/mesh_info.json";
    try {
        info.Write(record_file);
    } catch (const std::exception &e) {
        LOGF_ERR("Could not write to file " << record_file << ". Exception: " << e.what());
    }

}

void ComputeAndSaveMinWeightTriangulation(
    const ProgramOptionVariables &pov,
    const std::string &output_dir,
    const std::string &type,
    const trivis::mesh::TriMesh &mesh,
    const std::vector<std::vector<double>> &weights,
    double time_weights,
    std::optional<double> vis_radius,
    long random_seed,
    std::optional<int> iter_limit_opt,
    std::optional<double> improvement_ratio_limit_opt,
    std::optional<double> time_limit_opt,
    std::optional<int> simple_poly_max_size_opt,
    dr::MapDrawer *drawer = nullptr
) {
    auto rng = std::mt19937(random_seed);
    trivis::utils::SimpleClock clock;

    const auto *mesh_prev_ptr = &mesh;
    double orig_weight = mesh_optim::ComputeMeshWeight(*mesh_prev_ptr, weights);
    LOGF_INF("Orig weight: " << orig_weight << ".");

    double prev_weight = orig_weight;
    double new_weight = prev_weight;

    int i = 0;
    MeshWithInfo mesh_res;
    while (!iter_limit_opt || iter_limit_opt > i++) {
        if (time_limit_opt && clock.TimeInSeconds() > time_limit_opt) {
            break;
        }
        auto mesh_new = mesh_optim::MinWeightTriangulationIteration(*mesh_prev_ptr, weights, simple_poly_max_size_opt, rng);
        new_weight = mesh_optim::ComputeMeshWeight(mesh_new, weights);
        double improvement_ratio = (prev_weight - new_weight) / std::abs(prev_weight);
        // LOGF_INF("[" << i << "] New weight: " << new_weight << "; improvement ratio: " << improvement_ratio << ".");
        if (improvement_ratio_limit_opt && improvement_ratio < improvement_ratio_limit_opt) {
            break;
        }
        mesh_res = MeshWithInfo{type, std::move(mesh_new), time_weights, clock.TimeInSeconds(), random_seed, new_weight, vis_radius, i};
        if (!pov.save_only_last_iteration) {
            SaveMeshWithInfo(pov, mesh_res, output_dir + "_it-" + std::to_string(i), drawer);
        }
        // warning: must no longer access mesh_new (it was moved)!
        mesh_prev_ptr = &mesh_res.mesh;
        prev_weight = new_weight;
    }
    if (pov.save_only_last_iteration) {
        SaveMeshWithInfo(pov, mesh_res, output_dir, drawer);
    }

    double percentage_improvement = (orig_weight - new_weight) / std::abs(prev_weight);
    LOGF_INF("[DONE] New weight: " << new_weight << "; improvement ratio: " << percentage_improvement << ".");
}

bool SaveWeights(
    const std::string &file,
    const std::vector<std::vector<double>> &weights,
    double time_weights
) {
    std::fstream fs;
    fs.open(file, std::ios::out | std::ios::trunc);
    if (!fs.is_open()) {
        LOGF_ERR("Cannot open file " << file << " for writing.");
        return false;
    }
    int n = static_cast<int>(weights.size());
    fs << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_weights << "\n";
    fs << n << "\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (weights[i][j] != std::numeric_limits<double>::max()) {
                fs << i << " " << j << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << weights[i][j] << "\n";
            }
        }
    }
    fs.close();
    return true;
}

bool LoadWeights(
    const std::string &file,
    std::vector<std::vector<double>> &weights,
    double &time_weights
) {
    std::fstream fs;
    fs.open(file, std::ios::in);
    if (!fs.is_open()) {
        LOGF_ERR("Cannot open file " << file << " for reading.");
        return false;
    }
    try {
        int n;
        fs >> time_weights >> n;
        weights = std::vector<std::vector<double>>(n, std::vector<double>(n, std::numeric_limits<double>::max()));
        while (!fs.eof()) {
            int i, j;
            double w;
            fs >> i >> j >> w;
            weights[i][j] = w;
        }
        fs.close();
        return true;
    } catch (const std::exception &e) {
        LOGF_ERR("Error while reading file " << file << ": " << e.what());
        fs.close();
        return false;
    }
}

void ConstructAllTriangulations(
    const ProgramOptionVariables &pov,
    const trivis::Trivis &vis,
    const std::string &output_dir
) {

    // Create a helper object for drawing.
    auto drawer = dr::MakeMapDrawer(vis.map());

    if (pov.draw) {
        // Draw the map.
        drawer.OpenPDF(output_dir + "/map.pdf");
        drawer.DrawMap();
        drawer.Close();
    }

    trivis::utils::SimpleClock clock;
    std::set<std::string> m_set(pov.methods.begin(), pov.methods.end());
    double time_weights = 0.0;
    std::vector<std::vector<double>> weights;
    long seed = pov.random_seed_method;

    const auto &samp_dist = pov.max_sample_dist;
    const auto &eps = pov.eps_dist_diff_collinear;
    const auto &long_edge_penalty = pov.weights_long_edge_penalty;
    bool optimize_weights = long_edge_penalty >= 0;
    std::optional<double> long_edge_thresh_opt = pov.weights_long_edge_threshold >= 0.0 ? std::make_optional(pov.weights_long_edge_threshold) : std::nullopt;
    std::optional<double> short_edge_ver_max_opt = pov.weights_shortest_edges_ver_max > 0 ? std::make_optional(pov.weights_shortest_edges_ver_max) : std::nullopt;

    bool w_load = pov.load_weights;
    bool w_save = pov.save_weights;
    bool w_precompute = pov.precompute_weights;

    dr::MapDrawer *dr_ptr = pov.draw ? &drawer : nullptr;
    std::optional<int> mwt_iter = pov.mwt_iter_limit > 0 ? std::make_optional(pov.mwt_iter_limit) : std::nullopt;
    std::optional<double> mwt_imp_ratio = pov.mwt_improvement_ratio_limit > 0.0 ? std::make_optional(pov.mwt_improvement_ratio_limit) : std::nullopt;
    std::optional<double> mwt_time = pov.mwt_time_limit > 0.0 ? std::make_optional(pov.mwt_time_limit) : std::nullopt;
    std::optional<int> mwt_max_size = pov.mwt_simple_polygon_max_size > 0.0 ? std::make_optional(pov.mwt_simple_polygon_max_size) : std::nullopt;

    std::string type;

    auto ComputeAndSaveTriangulation = [
        pov, output_dir, vis, seed, mwt_iter, mwt_imp_ratio, mwt_time, mwt_max_size, dr_ptr
    ](
        const std::string &name,
        const std::string &type,
        const std::vector<std::vector<double>> &weights,
        double time_weights,
        std::optional<double> vis_radius
    ) {
        LOGF_INF("Computing " << name << " mesh ...");
        const auto &mesh = vis.mesh();
        std::string output_dir_full = output_dir + "/" + name;
        double time_construction;
        if (type == "CDT") {
            time_weights = 0.0;
            trivis::utils::SimpleClock clock;
            auto cdt = mesh_optim::GenerateConstrainedDelaunayTriangulation(vis);
            time_construction = clock.TimeInSeconds();
            auto mwi = MeshWithInfo{type, std::move(cdt), time_weights, time_construction};
            SaveMeshWithInfo(pov, mwi, output_dir + "/" + name, dr_ptr);
        } else {
            trivis::utils::SimpleClock clock;
            ComputeAndSaveMinWeightTriangulation(
                pov, output_dir_full, type, mesh, weights, time_weights, vis_radius, seed, mwt_iter, mwt_imp_ratio, mwt_time, mwt_max_size, dr_ptr);
            time_construction = clock.TimeInSeconds();
        }
        LOGF_INF("FINISHED | time weights: " << time_weights << ", time construction " << time_construction);
    };

    ///====================================///
    /// Constrained Delaunay Triangulation ///
    ///====================================///
    type = "CDT";
    if (Contains(m_set, type)) {
        const std::string &name = type;
        if (!w_precompute) {
            ComputeAndSaveTriangulation(name, type, weights, time_weights, std::nullopt);
        }
    }
    ///====================================///

    const std::string weights_file_template =
        pov.weights_dir + "/weights_" + pov.map_name + "_scale-" + std::to_string(pov.map_scale) +
            "_bucket-" + std::to_string(pov.bucket_size) + "-" + std::to_string(pov.mean_bucket_triangle_count_max) +
            "_eps-" + std::to_string(eps);

    auto GetWeights = [
        w_load, w_save, w_precompute
    ](
        const std::string &name,
        const std::string &weights_file,
        const std::function<double(std::vector<std::vector<double>> &)> &compute_weights_fun,
        std::vector<std::vector<double>> &weights,
        double &time_weights
    ) {
        bool weights_loaded = false;
        if (w_load && !w_precompute) {
            weights_loaded = LoadWeights(weights_file, weights, time_weights);
            if (!weights_loaded) {
                LOGF_WRN("Could not load weights from file " << weights_file << ".");
            }
        }
        if (!weights_loaded) {
            LOGF_INF("Computing " << name << " weights ...");
            time_weights = compute_weights_fun(weights);
            LOGF_INF("Done: " << time_weights << " s.");
            if (w_save || w_precompute) {
                bool ok = SaveWeights(weights_file, weights, time_weights);
                if (!ok) {
                    LOGF_WRN("Could not save weights to file " << weights_file << ".");
                }
            }
        }
    };

    ///==============================///
    /// Minimum Length Triangulation ///
    ///==============================///
    type = "MinLT";
    if (Contains(m_set, type)) {
        const std::string &name = type;
        std::string weights_file = weights_file_template + "_" + name + ".txt";
        auto ComputeWeights = [
            vis, eps
        ](std::vector<std::vector<double>> &weights) {
            trivis::utils::SimpleClock clock;
            weights = mesh_optim::ComputeWeightsEdgeLength(vis, eps);
            return clock.TimeInSeconds();
        };
        GetWeights(name, weights_file, ComputeWeights, weights, time_weights);
        if (!w_precompute) {
            ComputeAndSaveTriangulation(name, type, weights, time_weights, std::nullopt);
        }
    }
    ///==============================///
    /// Maximum Length Triangulation ///
    ///==============================///
    type = "MaxLT";
    if (Contains(m_set, type)) {
        const std::string &name = type;
        std::string weights_file = weights_file_template + "_" + name + ".txt";
        auto ComputeWeights = [
            vis, eps, m_set, time_weights
        ](std::vector<std::vector<double>> &weights) {
            if (Contains(m_set, "MinLT")) {
                // Reuse the weights from MinLT.
                trivis::utils::SimpleClock clock;
                mesh_optim::NegateWeights(weights);
                return time_weights + clock.TimeInSeconds();
            } else {
                trivis::utils::SimpleClock clock;
                weights = mesh_optim::ComputeWeightsEdgeLength(vis, eps);
                mesh_optim::NegateWeights(weights);
                return clock.TimeInSeconds();
            }
        };
        GetWeights(name, weights_file, ComputeWeights, weights, time_weights);
        if (!w_precompute) {
            ComputeAndSaveTriangulation(name, type, weights, time_weights, std::nullopt);
        }
    }
    ///==============================///

    auto GetWeightsFileVT = [
        weights_file_template,
        samp_dist,
        optimize_weights,
        long_edge_penalty,
        long_edge_thresh_opt,
        short_edge_ver_max_opt
    ](const std::string &name) {
        std::string ret = weights_file_template + "_" + name + "_sampdist-" + std::to_string(samp_dist);
        if (optimize_weights) {
            ret += "_penalty-" + std::to_string(long_edge_penalty);
            if (long_edge_thresh_opt) {
                ret += "-" + std::to_string(*long_edge_thresh_opt);
            }
            if (short_edge_ver_max_opt) {
                ret += "-ver-max-" + std::to_string(*short_edge_ver_max_opt);
            }
        }
        ret += ".txt";
        return ret;
    };

    ///==================================///
    /// Minimum Visibility Triangulation ///
    ///==================================///
    type = "MinVT";
    if (Contains(m_set, type)) {
        const std::string &name = type;
        std::string weights_file = GetWeightsFileVT(name);
        auto ComputeWeights = [
            vis, eps, optimize_weights, samp_dist, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt
        ](std::vector<std::vector<double>> &weights) {
            trivis::utils::SimpleClock clock;
            if (optimize_weights) {
                weights = mesh_optim::ComputeWeightsEdgeVisibilityOptimized(vis, samp_dist, eps, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt);
            } else {
                weights = mesh_optim::ComputeWeightsEdgeVisibility(vis, samp_dist, eps);
            }
            return clock.TimeInSeconds();
        };
        GetWeights(name, weights_file, ComputeWeights, weights, time_weights);
        if (!w_precompute) {
            ComputeAndSaveTriangulation(name, type, weights, time_weights, std::nullopt);
        }
    }
    ///==================================///
    /// Maximum Visibility Triangulation ///
    ///==================================///
    type = "MaxVT";
    if (Contains(m_set, "MaxVT")) {
        const std::string &name = type;
        std::string weights_file = GetWeightsFileVT(name);
        auto ComputeWeights = [
            vis, eps, optimize_weights, samp_dist, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt, m_set, time_weights
        ](std::vector<std::vector<double>> &weights) {
            if (Contains(m_set, "MinVT")) {
                // Reuse the weights from MinVT.
                trivis::utils::SimpleClock clock;
                mesh_optim::NegateWeights(weights);
                return clock.TimeInSeconds() + time_weights;
            } else {
                trivis::utils::SimpleClock clock;
                if (optimize_weights) {
                    weights = mesh_optim::ComputeWeightsEdgeVisibilityOptimized(vis, samp_dist, eps, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt);
                } else {
                    weights = mesh_optim::ComputeWeightsEdgeVisibility(vis, samp_dist, eps);
                }
                mesh_optim::NegateWeights(weights);
                return clock.TimeInSeconds();
            }
        };
        GetWeights(name, weights_file, ComputeWeights, weights, time_weights);
        if (!w_precompute) {
            ComputeAndSaveTriangulation(name, type, weights, time_weights, std::nullopt);
        }
    }
    ///==================================///

    for (double vis_radius: pov.vis_radii_mesh) {

        if (vis_radius <= 0.0) {
            LOGF_WRN("Ignoring vis. radius smaller or equal to zero.");
            continue;
        }

        ///====================================///
        /// Minimum d-Visibility Triangulation ///
        ///====================================///
        type = "d-MinVT";
        if (Contains(m_set, type)) {
            std::string name = std::to_string(vis_radius) + "-MinVT";
            std::string weights_file = GetWeightsFileVT(name);
            auto ComputeWeights = [
                vis, eps, optimize_weights, samp_dist, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt, vis_radius
            ](std::vector<std::vector<double>> &weights) {
                trivis::utils::SimpleClock clock;
                if (optimize_weights) {
                    weights = mesh_optim::ComputeWeightsEdgeVisibilityOptimized(vis, samp_dist, eps, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt, vis_radius);
                } else {
                    weights = mesh_optim::ComputeWeightsEdgeVisibility(vis, samp_dist, eps, vis_radius);
                }
                return clock.TimeInSeconds();
            };
            GetWeights(name, weights_file, ComputeWeights, weights, time_weights);
            if (!w_precompute) {
                ComputeAndSaveTriangulation(name, type, weights, time_weights, vis_radius);
            }
        }
        ///====================================///
        /// Maximum d-Visibility Triangulation ///
        ///====================================///
        type = "d-MaxVT";
        if (Contains(m_set, "d-MaxVT")) {
            std::string name = std::to_string(vis_radius) + "-MaxVT";
            std::string weights_file = GetWeightsFileVT(name);
            auto ComputeWeights = [
                vis, eps, optimize_weights, samp_dist, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt, vis_radius, m_set, time_weights
            ](std::vector<std::vector<double>> &weights) {
                if (Contains(m_set, "d-MinVT")) {
                    // Reuse the weights from d-MinVT.
                    trivis::utils::SimpleClock clock;
                    mesh_optim::NegateWeights(weights);
                    return clock.TimeInSeconds() + time_weights;
                } else {
                    trivis::utils::SimpleClock clock;
                    if (optimize_weights) {
                        weights = mesh_optim::ComputeWeightsEdgeVisibilityOptimized(vis, samp_dist, eps, long_edge_penalty, long_edge_thresh_opt, short_edge_ver_max_opt, vis_radius);
                    } else {
                        weights = mesh_optim::ComputeWeightsEdgeVisibility(vis, samp_dist, eps, vis_radius);
                    }
                    mesh_optim::NegateWeights(weights);
                    return clock.TimeInSeconds();
                }
            };
            GetWeights(name, weights_file, ComputeWeights, weights, time_weights);
            if (!w_precompute) {
                ComputeAndSaveTriangulation(name, type, weights, time_weights, vis_radius);
            }
        }
        ///====================================///
    }
}

void PrintMapInfo(const ProgramOptionVariables &pov, const trivis::Trivis &vis) {
    double area = 0.0;
    for (int i = 0; i < vis.triangles().size(); ++i) {
        const auto &triangle = vis.triangles()[i];
        area += trivis::geom::Area(triangle);
    }
    std::cout << pov.map_name;
    std::cout << " & " << pov.map_name;
    std::cout << " & " << vis.mesh().vertices.size();
    std::cout << " & " << vis.map().holes().size();
    std::cout << " & " << std::fixed << std::setprecision(2) << vis.limits().x_max - vis.limits().x_min;
    std::cout << " & " << std::fixed << std::setprecision(2) << vis.limits().y_max - vis.limits().y_min;
    std::cout << " & " << std::fixed << std::setprecision(2) << area;
    std::cout << " & ";
    std::cout << "\\\\\n";
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

    // Create and initialize Trivis object.
    trivis::Trivis vis;
    {   // Load map from file and move it to Trivis (without copying).
        std::stringstream info;
        auto map = trivis_plus::data_loading::LoadPolyMap(pov.map_full_path, pov.map_scale, &info);
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
    if (pov.mean_bucket_triangle_count_max > 0.0) {
        vis.FillPointLocationBuckets(std::nullopt, pov.mean_bucket_triangle_count_max);
    } else {
        vis.FillPointLocationBuckets(pov.bucket_size);
    }
    vis.OptimizePointLocationBucketTriangles();

    if (pov.print_map_info) {
        PrintMapInfo(pov, vis);
        exit(EXIT_SUCCESS);
    }

    // Create the output directory.
    std::string output_dir = pov.output_dir + "/" + pov.map_name + "_scale-" + std::to_string(pov.map_scale);
    fs::create_directories(output_dir);

    if (pov.precompute_weights || pov.save_weights) {
        fs::create_directories(pov.weights_dir);
    }

    ConstructAllTriangulations(pov, vis, output_dir);

    if (pov.precompute_weights) {
        return EXIT_SUCCESS;
    }

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
