/**
 * File:   gen_points.cc
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
#include "trivis_plus/data_loading/load_mesh.h"
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

struct ProgramOptionVariables {
    std::string map_name;
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string mesh_dir = DEFAULT_MESH_DIR;
    std::string points_dir = DEFAULT_POINT_DIR;
    double bucket_size = 1.0;
    int max_points = 1000;
    int eps_close_point_exp_max = -1;
    int eps_close_point_exp_min = -15;
    int random_seed = 42;
    bool draw = false;
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
        ("random-seed", po::value(&pov.random_seed)->default_value(pov.random_seed), "Random seed.")
        ("draw", po::bool_switch(&pov.draw)->default_value(pov.draw), "Draw the mesh and the points.");
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
    return '0';
}

bool SaveAndDrawPoints(
    const trivis::geom::FPoints &points,
    const std::string &points_path_no_ext,
    const std::string &points_name,
    dr::MapDrawer &drawer,
    bool draw = true
) {
    trivis::utils::SimpleClock clock;
    LOGF_INF("Saving " << points.size() << " " << points_name << " points.");
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
    if (draw) {
        drawer.OpenPDF(points_path_no_ext + "_" + points_name + ".pdf");
        drawer.DrawMap();
        drawer.DrawPoints(points, 0.5, dr::kColorRed);
        drawer.Close();
    }
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
    return true;
}

int MainBody(const ProgramOptionVariables &pov) {

    trivis::utils::SimpleClock clock;
    std::stringstream info;

    LOGF_INF("Loading map from " << pov.map_full_path << ".");
    clock.Restart();
    auto map = trivis_plus::data_loading::LoadPolyMap(pov.map_full_path, std::nullopt, &info);
    if (!map) {
        LOGF_FTL("Error while loading map: '" << info.str() << "'.");
        return EXIT_FAILURE;
    }
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Initializing TřiVis.");
    clock.Restart();
    trivis::Trivis vis{std::move(map.value())}; // Default initialization.
    map = std::nullopt; // cannot use map anymore ! (it was moved)
    const auto &lim = vis.limits();
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
    LOGF_INF("Map limits: [ MIN: " << lim.x_min << ", " << lim.y_min << " | MAX: " << lim.x_max << ", " << lim.y_max << " ].");

    LOGF_INF("Saving mesh.");
    clock.Restart();
    std::string mesh_path_no_ext = pov.mesh_dir + "/" + pov.map_name + "_cdt";
    trivis_plus::data_loading::SaveTriMesh(vis.mesh(), mesh_path_no_ext + ".txt");
    auto drawer = dr::MakeMapDrawer(vis.map());
    if (pov.draw) {
        drawer.OpenPDF(mesh_path_no_ext + ".pdf");
        drawer.DrawMap();
        drawer.DrawPolygons(vis.triangles(), 0.1, dr::kColorRed);
        drawer.Close();
    }
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Generating points.");
    clock.Restart();
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

    {   // Generate points on map vertices (Ver).
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
        if (!SaveAndDrawPoints(points, point_path_no_ext, "Ver", drawer, pov.draw)) {
            return EXIT_FAILURE;
        }
    }

    {   // Generate points close to map vertices (NearV).
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
        if (!SaveAndDrawPoints(points, point_path_no_ext, "NearV", drawer, pov.draw)) {
            return EXIT_FAILURE;
        }
    }

    {   // Generate points on map edge midpoints (Mid).
        auto rng = std::mt19937(pov.random_seed);
        points.clear();
        points.reserve(vis.mesh().edges.size());
        for (const auto &edge: vis.mesh().edges) {
            auto midpoint = (vis.mesh().vertices[edge.vertices[0]].point + vis.mesh().vertices[edge.vertices[1]].point) / 2.0;
            points.push_back(std::move(midpoint));
        }
        std::shuffle(points.begin(), points.end(), rng);
        if (points.size() > pov.max_points) {
            points.resize(pov.max_points);
        }
        if (!SaveAndDrawPoints(points, point_path_no_ext, "Mid", drawer, pov.draw)) {
            return EXIT_FAILURE;
        }
    }

    {   // Generate points close to map edge midpoints (NearM).
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
        if (!SaveAndDrawPoints(points, point_path_no_ext, "NearM", drawer, pov.draw)) {
            return EXIT_FAILURE;
        }
    }

    {   // Generate points inside the bounding box of the map (BB).
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
        if (!SaveAndDrawPoints(points, point_path_no_ext, "BB", drawer, pov.draw)) {
            return EXIT_FAILURE;
        }
    }

    {   // Generate points inside map (In).
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
        if (!SaveAndDrawPoints(points, point_path_no_ext, "In", drawer, pov.draw)) {
            return EXIT_FAILURE;
        }
    }
    LOGF_INF("DONE (Generating points). It took " << clock.TimeInSeconds() << " seconds.");

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
