/**
 * File:   generate_instances.cc
 *
 * Date:   10.01.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "trivis/trivis.h"
#include "trivis/utils/random_points.h"
#include "trivis/utils/simple_clock.h"

#include "trivis_plus/data_loading/load_map.h"
#include "trivis_plus/drawing/fancy_drawing.h"
#include "trivis_plus/utils/log.h"

#ifndef DEFAULT_MAP_DIR
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_MESH_DIR
#define DEFAULT_MESH_DIR "."
#endif

#ifndef DEFAULT_POINT_DIR
#define DEFAULT_POINT_DIR "."
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;
namespace dr = trivis_plus::drawing;

/**
 * All program option variables and their default values should be defined here.
 * For each variable, there should be an option added in AddProgramOptions.
 */
struct ProgramOptionVariables {
    std::string map_name;
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string mesh_dir = DEFAULT_MESH_DIR;
    std::string points_dir = DEFAULT_POINT_DIR;
    double bucket_size = 1.0;
    int max_points = 10000;
    int eps_close_point_exp_max = -1;
    int eps_close_point_exp_min = -15;
    int random_seed = 42;
};

void AddProgramOptions(
    po::options_description &options_description,
    ProgramOptionVariables &pov
) {
    options_description.add_options()
        ("help,h", "Produce this help message. \n (*) Overwrites options: all.")
        ("map-name", po::value(&pov.map_name)->default_value(pov.map_name), "Map name.")
        ("map", po::value(&pov.map_full_path)->default_value(pov.map_full_path), "Full path to the map file. \n (*) Overwrites options: map-name, map-ext, map-dir.")
        ("map-ext", po::value(&pov.map_extension)->default_value(pov.map_extension), "Map file extension.")
        ("map-dir", po::value(&pov.map_dir)->default_value(pov.map_dir), "Map file directory.")
        ("mesh-dir", po::value(&pov.mesh_dir)->default_value(pov.mesh_dir), "Mesh file directory.")
        ("points-dir", po::value(&pov.points_dir)->default_value(pov.points_dir), "Point file directory.")
        ("bucket-size", po::value(&pov.bucket_size)->default_value(pov.bucket_size), "Bucket size.")
        ("max-points", po::value(&pov.max_points)->default_value(pov.max_points), "Number of points.")
        ("eps-close-point-exp-max", po::value(&pov.eps_close_point_exp_max)->default_value(pov.eps_close_point_exp_max), "Max exponent of epsilon for close points.")
        ("eps-close-point-exp-min", po::value(&pov.eps_close_point_exp_min)->default_value(pov.eps_close_point_exp_min), "Min exponent of epsilon for close points.")
        ("random-seed", po::value(&pov.random_seed)->default_value(pov.random_seed), "Random seed.");
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
    return '0';
}

bool SaveAndDrawPoints(
    const trivis::geom::FPoints &points,
    const std::string &points_path_no_ext,
    const std::string &points_name,
    dr::MapDrawer &drawer
) {
    trivis::utils::SimpleClock clock;
    LOGF_INF(">> Saving " << points.size() << " " << points_name << " points and drawing.");
    clock.Restart();
    std::fstream fs;
    fs.open(points_path_no_ext + "_" + points_name + ".txt", std::ios::out | std::ios::trunc);
    if (!fs.is_open()) {
        LOGF_FTL("Cannot open file " << points_path_no_ext << ".txt for writing.");
        return false;
    }
    fs << points.size() << "\n";
    for (int i = 0; i < points.size(); ++i) {
        const auto &point = points[i];
        fs << i << " "
           << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << point.x << " "
           << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << point.y << "\n";
    }
    fs.close();
    drawer.OpenPDF(points_path_no_ext + "_" + points_name + ".pdf");
    drawer.DrawMap();
    drawer.DrawPoints(points, 0.25, dr::kColorRed);
    drawer.Close();
    LOGF_INF("<< DONE. It took " << clock.TimeInSeconds() << " seconds.");
    return true;
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

    trivis::utils::SimpleClock clock;

    LOGF_INF(">> Initializing TriVis.");
    clock.Restart();
    // Create and initialize the TriVis object.
    trivis::Trivis vis;
    {   // Load map from file and move it to TriVis (without copying).
        trivis::geom::PolyMap map;
        std::string load_msg = trivis_plus::data_loading::LoadPolyMapSafely(pov.map_full_path, map);
        if (load_msg != "ok") {
            LOGF_FTL("Error while loading map. " << load_msg);
            return EXIT_FAILURE;
        }
        vis.SetMap(std::move(map)); // Set the map to TriVis instance.
    }
    vis.ConstructMeshCDT();
    LOGF_INF("<< DONE. It took " << clock.TimeInSeconds() << " seconds.");

    const auto &lim = vis.limits();
    LOGF_INF("Map limits [ MIN: " << lim.x_min << ", " << lim.y_min << " | MAX: " << lim.x_max << ", " << lim.y_max << " ].");

    LOGF_INF(">> Saving mesh and drawing.");
    clock.Restart();
    std::string mesh_path_no_ext = pov.mesh_dir + "/" + pov.map_name + "_cdt";
    trivis::Trivis::SaveMesh(vis.mesh(), mesh_path_no_ext + ".txt");
    auto drawer = dr::MakeMapDrawer(vis.map());
    drawer.OpenPDF(mesh_path_no_ext + ".pdf");
    drawer.DrawMap();
    drawer.DrawPolygons(vis.triangles(), 0.05, dr::kColorRed);
    drawer.Close();
    LOGF_INF("<< DONE. It took " << clock.TimeInSeconds() << " seconds.");

    int n_close_points = pov.eps_close_point_exp_max - pov.eps_close_point_exp_min + 1;
    std::vector<std::normal_distribution<double>> normal_dists;
    normal_dists.reserve(n_close_points);
    for (int i = 0; i < n_close_points; ++i) {
        normal_dists.emplace_back(0.0, std::pow(10.0, pov.eps_close_point_exp_max - i));
    }
    std::uniform_real_distribution<double> uniform_dist_x(lim.x_min, lim.x_max);
    std::uniform_real_distribution<double> uniform_dist_y(lim.y_min, lim.y_max);

    trivis::geom::FPoints points;
    auto point_path_no_ext = pov.points_dir + "/" + pov.map_name + "_points";

    {   // Get node points.
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(vis.mesh().vertices.size());
        for (const auto &node: vis.mesh().vertices) {
            points.push_back(node.point);
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "on_nodes", drawer)) {
            return EXIT_FAILURE;
        }
    }

    {   // Get close-to-node points.
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(vis.mesh().vertices.size() * n_close_points);
        for (const auto &node: vis.mesh().vertices) {
            for (int i = 0; i < n_close_points; ++i) {
                points.emplace_back(node.point.x + normal_dists[i](rng), node.point.y + normal_dists[i](rng));
            }
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "close_to_nodes", drawer)) {
            return EXIT_FAILURE;
        }
    }

    {   // Get edge mid-points.
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(vis.mesh().edges.size());
        for (const auto &edge: vis.mesh().edges) {
            points.emplace_back((vis.mesh().vertices[edge.vertices[0]].point + vis.mesh().vertices[edge.vertices[1]].point) / 2.0);
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "edge_mid_points", drawer)) {
            return EXIT_FAILURE;
        }
    }

    {   // Get close-to-edge-mid-points.
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(vis.mesh().edges.size() * n_close_points);
        for (const auto &edge: vis.mesh().edges) {
            for (int i = 0; i < n_close_points; ++i) {
                auto midpoint = (vis.mesh().vertices[edge.vertices[0]].point + vis.mesh().vertices[edge.vertices[1]].point) / 2.0;
                points.emplace_back(midpoint.x + normal_dists[i](rng), midpoint.y + normal_dists[i](rng));
            }
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "close_to_edge_mid_points", drawer)) {
            return EXIT_FAILURE;
        }
    }

    {   // Get random points.
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(pov.max_points);
        for (int i = 0; i < pov.max_points; ++i) {
            points.emplace_back(uniform_dist_x(rng), uniform_dist_y(rng));
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "random", drawer)) {
            return EXIT_FAILURE;
        }
    }

    {   // Get random interior points.
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(pov.max_points);
        for (int i = 0; i < pov.max_points; ++i) {
            auto point = trivis::utils::UniformRandomPointInRandomTriangle(vis.triangles(), rng);
            points.push_back(point);
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "random_interior", drawer)) {
            return EXIT_FAILURE;
        }
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
