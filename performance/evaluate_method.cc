/**
 * File:   evaluate_method.cc
 *
 * Date:   22.01.2024
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
#include "trivis/utils/clipper_geom.h"

#include "trivis_plus/data_loading/load_map.h"
#include "trivis_plus/utils/log.h"
#include "trivis_plus/drawing/fancy_drawing.h"

#ifndef DEFAULT_MAP_DIR
#define DEFAULT_MAP_DIR "."
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
    std::string map_type;
    std::string map_name;
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string points_type;
    std::string points_name;
    std::string points_extension = ".txt";
    std::string points_dir = DEFAULT_POINT_DIR;
    std::string points_full_path;
    std::string algorithm_name;
    std::string algorithm_data_file;
    std::string algorithm_gt_name;
    std::string algorithm_gt_data_file;
    std::string out_file;
    bool out_file_append = false;
};

template<typename T>
inline auto po_value(T &variable) {
    return po::value(&variable)->default_value(variable);
}

template<typename T>
inline auto po_values(std::vector<T> &variables) {
    assert(variables.empty());
    return po::value(&variables)->multitoken()->default_value(variables, "[]");
}

void AddProgramOptions(
    po::options_description &options_description,
    ProgramOptionVariables &pov
) {
    options_description.add_options()
        ("help,h", "Produce this help message. \n (*) Overwrites options: all.")
        ("map-type", po_value(pov.map_type), "Map type.")
        ("map-name", po_value(pov.map_name), "Map name.")
        ("map-extension", po_value(pov.map_extension), "Map extension.")
        ("map-dir", po_value(pov.map_dir), "Map directory.")
        ("map-full-path", po_value(pov.map_full_path), "Full path to the map file. \n (*) Overwrites options: map-name, map-ext, map-dir.")
        ("points-type", po_value(pov.points_type), "Points type.")
        ("points-name", po_value(pov.points_name), "Points name.")
        ("points-extension", po_value(pov.points_extension), "Points extension.")
        ("points-dir", po_value(pov.points_dir), "Points directory.")
        ("points-full-path", po_value(pov.points_full_path), "Full path to the points file. \n (*) Overwrites options: points-name, points-ext, points-dir.")
        ("algorithm-name", po_value(pov.algorithm_name), "Algorithm name.")
        ("algorithm-data-file", po_value(pov.algorithm_data_file), "Algorithm data file.")
        ("algorithm-gt-name", po_value(pov.algorithm_gt_name), "Algorithm ground truth name.")
        ("algorithm-gt-data-file", po_value(pov.algorithm_gt_data_file), "Algorithm ground truth data file.")
        ("out-file", po_value(pov.out_file), "Output file.")
        ("out-file-append", po_value(pov.out_file_append), "Append to output file.");
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

struct LoadPolygonResult {
    bool crashed = false;
    bool timeout = false;
    double time_locate = std::numeric_limits<double>::quiet_NaN();
    double time_query = std::numeric_limits<double>::quiet_NaN();
    double time_intersect = std::numeric_limits<double>::quiet_NaN();
    double time_query_intersect = std::numeric_limits<double>::quiet_NaN();
    double time_total = std::numeric_limits<double>::quiet_NaN();
    bool snapped = false;
    bool found = false;
    trivis::geom::FPolygon polygon;
};

std::optional<LoadPolygonResult> LoadPolygon(
    int query_id,
    const trivis::geom::FPoint &query,
    std::ifstream &ifs
) {
    LoadPolygonResult result;
    std::string line;
    std::string token;
    std::istringstream iss;
    std::getline(ifs, line);
    iss = std::istringstream(line);
    int id_check;
    iss >> id_check;
    if (id_check != query_id) {
        LOGF_FTL("Point ids in algorithm data files are not consistent: " << id_check << " != " << query_id << " on line '" << line << "'.");
        return std::nullopt;
    }
    iss >> token;
    if (token == "crashed") {
        result.crashed = true;
        return result;
    } else if (token == "timeout") {
        result.timeout = true;
        return result;
    } else if (token != "times") {
        LOGF_FTL("Unexpected token '" << token << "' on line '" << line << "'.");
        return std::nullopt;
    }
    int n_times = 0;
    iss >> n_times;
    std::vector<double> times(n_times, 0.0);
    for (int k = 0; k < n_times; ++k) {
        iss >> times[k];
    }
    if (n_times == 1) {
        result.time_total = times[0];
    } else if (n_times == 2) {
        result.time_locate = times[0];
        result.time_query_intersect = times[1];
        result.time_total = result.time_locate + result.time_query_intersect;
    } else if (n_times == 3) {
        result.time_locate = times[0];
        result.time_query = times[1];
        result.time_intersect = times[2];
        result.time_query_intersect = result.time_query + result.time_intersect;
        result.time_total = result.time_locate + result.time_query_intersect;
    } else if (n_times >= 4) {
        LOGF_FTL("Unexpected number of times '" << n_times << "' on line '" << line << "'.");
        return std::nullopt;
    }
    iss >> token;
    std::optional<trivis::geom::FPoint> snap_point;
    if (token == "snapped") {
        double x, y;
        iss >> x >> y;
        snap_point = trivis::geom::FPoint(x, y);
        iss >> token;
    }
    if (snap_point && *snap_point != query) {
        result.snapped = true;
    }
    if (token != "polygon") {
        LOGF_FTL("Unexpected token (polygon) '" << token << "' on line '" << line << "'.");
        return std::nullopt;
    }
    int n_polygon;
    iss >> n_polygon;
    if (n_polygon > 0) {
        result.found = true;
    }
    result.polygon.resize(n_polygon);
    for (int k = 0; k < n_polygon; ++k) {
        iss >> result.polygon[k].x >> result.polygon[k].y;
    }
    return result;
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
    std::string line;

    trivis::geom::PolyMap map;
    std::string load_msg = trivis_plus::data_loading::LoadPolyMapSafely(pov.map_full_path, map);
    if (load_msg != "ok") {
        LOGF_FTL("Error while loading map. " << load_msg);
        return EXIT_FAILURE;
    }
    trivis::geom::FLimits lim = map.limits();
    Clipper2Lib::Paths64 map_clipper = trivis::utils::ToClipper(map);
    double map_area_clipper = Clipper2Lib::Area(map_clipper);

    trivis::Trivis vis;
    vis.SetMap(map);
    vis.ConstructMeshCDT();
    vis.FillPointLocationBuckets(1.0);

    std::ifstream ifs_points(pov.points_full_path);
    if (!ifs_points.is_open()) {
        LOGF_FTL("Cannot open file " << pov.points_full_path << " for reading.");
        return EXIT_FAILURE;
    }

    // Skip the first line in the points file.
    std::getline(ifs_points, line);

    std::ifstream ifs_algorithm_data(pov.algorithm_data_file);
    if (!ifs_algorithm_data.is_open()) {
        LOGF_FTL("Cannot open file " << pov.algorithm_data_file << " for reading.");
        return EXIT_FAILURE;
    }

    std::ifstream ifs_algorithm_gt_data(pov.algorithm_gt_data_file);
    if (!ifs_algorithm_gt_data.is_open()) {
        LOGF_FTL("Cannot open file " << pov.algorithm_gt_data_file << " for reading.");
        return EXIT_FAILURE;
    }

    // Skip the first three lines in the ground truth algorithm data file.
    std::getline(ifs_algorithm_gt_data, line);
    std::getline(ifs_algorithm_gt_data, line);
    std::getline(ifs_algorithm_gt_data, line);

    std::ofstream ofs(pov.out_file, pov.out_file_append ? std::ios_base::app : std::ios_base::out);
    if (!ofs.is_open()) {
        LOGF_FTL("Cannot open file " << pov.out_file << " for writing.");
        return EXIT_FAILURE;
    }

    std::string token;
    std::istringstream iss;
    std::getline(ifs_algorithm_data, line);
    iss = std::istringstream(line);
    double time_init;
    iss >> token >> time_init;
    if (token != "time_init") {
        LOGF_FTL("Unexpected token " << token << ".");
        return EXIT_FAILURE;
    }
    std::getline(ifs_algorithm_data, line);
    iss = std::istringstream(line);
    double time_preprocess;
    iss >> token >> time_preprocess;
    if (token != "time_preprocess") {
        LOGF_FTL("Unexpected token " << token << ".");
        return EXIT_FAILURE;
    }
    std::getline(ifs_algorithm_data, line);
    iss = std::istringstream(line);
    int n_points;
    iss >> n_points;

    if (!pov.out_file_append) {
        ofs << "map_type";
        ofs << ",map_name";
        ofs << ",map_area";
        ofs << ",points";
        ofs << ",algorithm";
        ofs << ",time_init";
        ofs << ",time_preprocess";
        ofs << ",query_id";
        ofs << ",node_id";
        ofs << ",weakly_simple";
        ofs << ",crashed";
        ofs << ",timeout";
        ofs << ",time_locate";
        ofs << ",time_query";
        ofs << ",time_intersect";
        ofs << ",time_query_intersect";
        ofs << ",time_total";
        ofs << ",snapped";
        ofs << ",found";
        ofs << ",area";
        ofs << ",gt_algorithm";
        ofs << ",gt_available";
        ofs << ",gt_found";
        ofs << ",gt_area";
        ofs << ",gt_diff_area";
        ofs << ",diff_gt_area";
        ofs << "\n";
    }

    std::ostringstream oss;
    oss << pov.map_type;
    oss << "," << pov.map_name;
    oss << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << map_area_clipper;
    oss << "," << pov.points_type;
    oss << "," << pov.algorithm_name;
    oss << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_init;
    oss << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_preprocess;
    std::string common_prefix = oss.str();

    for (int i = 0; i < n_points; ++i) {
        // Load and locate the query point.
        int query_id;
        double x, y;
        std::getline(ifs_points, line);
        iss = std::istringstream(line);
        iss >> query_id >> x >> y;
        trivis::geom::FPoint query(x, y);
        int node_id = -1;
        bool weakly_simple = false;
        auto query_locate = vis.LocatePoint(query, 0.0);
        if (query_locate && !query_locate->snap_to_nodes.empty()) {
            node_id = query_locate->snap_to_nodes.front();
            weakly_simple = query_locate->snap_to_nodes.size() > 1;
        }
        // Load the algorithm data.
        auto algorithm_result = LoadPolygon(query_id, query, ifs_algorithm_data);
        if (!algorithm_result) {
            LOGF_FTL("Error while loading algorithm data.");
            return EXIT_FAILURE;
        }
        auto gt_algorithm_result = LoadPolygon(query_id, query, ifs_algorithm_gt_data);
        if (!gt_algorithm_result) {
            LOGF_FTL("Error while loading algorithm ground truth data.");
            return EXIT_FAILURE;
        }
        // Get the outputs.
        bool crashed = algorithm_result->crashed;
        bool timeout = algorithm_result->timeout;
        double time_locate = algorithm_result->time_locate;
        double time_query = algorithm_result->time_query;
        double time_intersect = algorithm_result->time_intersect;
        double time_query_intersect = algorithm_result->time_query_intersect;
        double time_total = algorithm_result->time_total;
        bool snapped = algorithm_result->snapped;
        bool found = algorithm_result->found;
        trivis::geom::FPolygon polygon = algorithm_result->polygon;
        Clipper2Lib::Path64 polygon_clipper = trivis::utils::ToClipper(polygon, lim);
        double area = Clipper2Lib::Area(polygon_clipper);
        std::string gt_algorithm = pov.algorithm_gt_name;
        bool gt_available = !gt_algorithm_result->crashed && !gt_algorithm_result->timeout;
        bool gt_found = gt_algorithm_result->found;
        trivis::geom::FPolygon gt_polygon = gt_algorithm_result->polygon;
        Clipper2Lib::Path64 gt_polygon_clipper = trivis::utils::ToClipper(gt_polygon, lim);
        double gt_area = Clipper2Lib::Area(gt_polygon_clipper);
        auto gt_diff = Clipper2Lib::Difference({gt_polygon_clipper}, {polygon_clipper}, Clipper2Lib::FillRule::NonZero);
        double gt_diff_area = Clipper2Lib::Area(gt_diff);
        auto diff_gt = Clipper2Lib::Difference({polygon_clipper}, {gt_polygon_clipper}, Clipper2Lib::FillRule::NonZero);
        double diff_gt_area = Clipper2Lib::Area(diff_gt);
        // Write the outputs.
        ofs << common_prefix;
        ofs << "," << query_id;
        ofs << "," << node_id;
        ofs << "," << weakly_simple;
        ofs << "," << crashed;
        ofs << "," << timeout;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_locate;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_query;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_intersect;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_query_intersect;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_total;
        ofs << "," << snapped;
        ofs << "," << found;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << area;
        ofs << "," << gt_algorithm;
        ofs << "," << gt_available;
        ofs << "," << gt_found;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << gt_area;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << gt_diff_area;
        ofs << "," << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << diff_gt_area;
        ofs << "\n";
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