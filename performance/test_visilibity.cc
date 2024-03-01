/**
 * File:   test_visilibity.cc
 *
 * Date:   11.01.2024
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

#include "visilibity/visilibity.hpp"

#ifndef DEFAULT_MAP_DIR
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_POINT_DIR
#define DEFAULT_POINT_DIR "."
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;

/**
 * All program option variables and their default values should be defined here.
 * For each variable, there should be an option added in AddProgramOptions.
 */
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
    int max_points = std::numeric_limits<int>::max();
    double epsilon = 0.0;
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
        ("max-points", po_value(pov.max_points), "Maximum number of points to test.")
        ("epsilon", po_value(pov.epsilon), "Epsilon.");
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

VisiLibity::Polygon ToVisiLibity(const trivis::geom::FPolygon &polygon) {
    VisiLibity::Polygon visi_polygon;
    for (const auto &point: polygon) {
        visi_polygon.push_back(VisiLibity::Point(point.x, point.y));
    }
    return visi_polygon;
}

std::vector<VisiLibity::Polygon> ToVisiLibity(const trivis::geom::PolyMap &map) {
    std::vector<VisiLibity::Polygon> visi_polygons;
    visi_polygons.push_back(ToVisiLibity(map.border()));
    for (const auto &hole: map.holes()) {
        visi_polygons.push_back(ToVisiLibity(hole));
    }
    return visi_polygons;
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

    if (pov.epsilon < 0.0) {
        LOGF_FTL("Epsilon must be non-negative.");
        return EXIT_FAILURE;
    }

    trivis::utils::SimpleClock clock;

    trivis::geom::PolyMap map;
    std::string load_msg = trivis_plus::data_loading::LoadPolyMapSafely(pov.map_full_path, map);
    if (load_msg != "ok") {
        LOGF_FTL("Error while loading map. " << load_msg);
        return EXIT_FAILURE;
    }
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

    auto visi_polygons = ToVisiLibity(map);

    clock.Restart();
    auto visi_env = VisiLibity::Environment(visi_polygons);
    double time_init = clock.TimeInSeconds();

    if (!pov.resume_after_crash) {
        ofs << "time_init" << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_init << "\n";
    }
    if (pov.check_map_validity) {
        if (!visi_env.is_valid(pov.epsilon)) {
            LOGF_WRN("Visi environment is not valid.");
        }
    }

    double time_preprocess = 0.0;

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
        auto query = VisiLibity::Point(x, y);
        bool is_inside = query.in(visi_env, pov.epsilon);
        double time_is_inside = clock.TimeInSeconds();
        auto query_old = query;
        clock.Restart();
        query.snap_to_boundary_of(visi_env, pov.epsilon);
        query.snap_to_vertices_of(visi_env, pov.epsilon);
        auto visi_polygon = is_inside ? std::make_optional<VisiLibity::Visibility_Polygon>(VisiLibity::Visibility_Polygon(query, visi_env, pov.epsilon)) : std::nullopt;
        double time_vis_poly = clock.TimeInSeconds();
        ofs << id << " times " << 2
            << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_is_inside
            << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << time_vis_poly;
        if (query != query_old) {
            ofs << " snapped"
                << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << query.x()
                << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << query.y();
        }
        ofs << " polygon ";
        if (visi_polygon) {
            ofs << visi_polygon->n();
            for (int j = 0; j < visi_polygon->n(); ++j) {
                ofs << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << (*visi_polygon)[j].x()
                    << " " << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10) << (*visi_polygon)[j].y();
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
    return MainBody(pov);
}