/**
 * File:   draw_visibility_region.cc
 *
 * Date:   04.11.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */


#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "trivis/trivis.h"

#include "trivis_plus/data_loading/load_map.h"
#include "trivis_plus/drawing/fancy_drawing.h"
#include "trivis_plus/utils/log.h"

#ifndef DEFAULT_MAP_DIR // May be defined in CMakeLists.txt.
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_OUT_DIR // May be defined in CMakeLists.txt.
#define DEFAULT_OUT_DIR "."
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;

/**
 * All program option variables and their default values should be defined here.
 * For each variable, there should be an option added in AddProgramOptions.
 */
struct ProgramOptionVariables {
    std::string verbosity = "info";
    std::string map_name = "potholes";
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string output_dir = DEFAULT_OUT_DIR;
    double map_scale = 0.01;
    bool shoot_ray = false;
    double x = 1.0;
    double y = 15.0;
    double tx = 19.0;
    double ty = 3.0;
    double vis_radius = -1.0;
};

void AddProgramOptions(
    po::options_description &options_description,
    ProgramOptionVariables &pov
) {
    options_description.add_options()
        ("help,h", "Produce this help message. \n (*) Overwrites options: all.")
        ("verbosity",
         po::value(&pov.verbosity)->default_value(pov.verbosity),
         "Log messages: fatal, error, warning, info, debug, trace.")
        ("map-name",
         po::value(&pov.map_name)->default_value(pov.map_name),
         "Map name.")
        ("map",
         po::value(&pov.map_full_path)->default_value(pov.map_full_path),
         "Full path to the map file. \n (*) Overwrites options: map-name, map-ext, map-dir.")
        ("map-ext",
         po::value(&pov.map_extension)->default_value(pov.map_extension),
         "Map file extension.")
        ("map-dir",
         po::value(&pov.map_dir)->default_value(pov.map_dir),
         "Map file directory.")
        ("out-dir",
         po::value(&pov.output_dir)->default_value(pov.output_dir),
         "Output directory.")
        ("map-scale",
         po::value(&pov.map_scale)->default_value(pov.map_scale),
         "Map coordinates are scaled by this factor when loading the map.")
        ("shoot-ray",
         po::bool_switch(&pov.shoot_ray)->default_value(pov.shoot_ray),
         "Shoot a ray from the source in the direction of the target.")
        ("x",
         po::value(&pov.x)->default_value(pov.x),
         "X coordinate of the source.")
        ("y",
         po::value(&pov.y)->default_value(pov.y),
         "Y coordinate of the source.")
        ("tx",
         po::value(&pov.tx)->default_value(pov.tx),
         "X coordinate of the target.")
        ("ty",
         po::value(&pov.ty)->default_value(pov.ty),
         "Y coordinate of the target.")
        ("vis-radius",
         po::value(&pov.vis_radius)->default_value(pov.vis_radius),
         "Limited visibility radius (-1 ~ infinite).");
}

/**
 * Get the severity level from the verbosity string.
 * @param pov Program option variables.
 * @return Severity level.
 */
trivis_plus::utils::severity_level GetSeverity(const ProgramOptionVariables &pov) {
    if (pov.verbosity == "trace") return trivis_plus::utils::severity_level::trace;
    if (pov.verbosity == "debug") return trivis_plus::utils::severity_level::debug;
    if (pov.verbosity == "info") return trivis_plus::utils::severity_level::info;
    if (pov.verbosity == "warning") return trivis_plus::utils::severity_level::warning;
    if (pov.verbosity == "error") return trivis_plus::utils::severity_level::error;
    if (pov.verbosity == "fatal") return trivis_plus::utils::severity_level::fatal;
    LOGF_WRN("Unknown verbosity level '" << pov.verbosity << "'. Using 'info' instead.");
    return trivis_plus::utils::severity_level::info;
}

/**
 * Parse the program options and initialize logging.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @param pov Program option variables.
 * @return Character 'e' if an exception occurred, 'h' if help option, '0' else.
 */
char ParseProgramOptions(
    int argc,
    const char *const *argv,
    ProgramOptionVariables &pov
) {
    po::variables_map vm;
    po::options_description command_line_options;
    po::options_description options_description("General options");
    AddProgramOptions(options_description, pov);
    try {
        // Parse the command line arguments.
        command_line_options.add(options_description);
        po::store(po::parse_command_line(argc, argv, command_line_options), vm);
        if (vm.count("help")) {
            // If help option, print the options and return 'h'.
            command_line_options.print(std::cout, 80);
            return 'h';
        }
        po::notify(vm);
    } catch (const std::exception &e) {
        // If exception, log it and return 'e'.
        trivis_plus::utils::InitLogging(trivis_plus::utils::severity_level::info);
        LOGF_FTL("Error in parsing arguments: " << e.what() << ".");
        return 'e';
    }
    trivis_plus::utils::InitLogging(GetSeverity(pov));
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

    LOGF_INF("Running the two-point visibility example.");

    trivis::utils::SimpleClock clock;
    std::stringstream info;

    LOGF_INF("Loading map from " << pov.map_full_path << ".");
    clock.Restart();
    auto map = trivis_plus::data_loading::LoadPolyMap(pov.map_full_path, pov.map_scale, &info);
    if (!map) {
        LOGF_FTL("Error while loading map: '" << info.str() << "'.");
        return EXIT_FAILURE;
    }
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Preprocessing the map.");
    clock.Restart();
    // Warning: order of the following matters!
    map->ShiftToOrigin(); // Subtracts min X and Y coordinates from all points.
    map->RemoveDuplicatePoints(); // Removes all consecutive identical points.
    map->SimplifyWeaklySelfIntersectingPolygons(); // Splits all weakly self-intersecting polygons to multiple (touching) simple polygons.
    map->RemoveCollinearPoints(); // Removes all consecutive collinear points.
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Initializing TřiVis.");
    clock.Restart();
    trivis::Trivis vis{std::move(map.value())}; // Default initialization.
    map = std::nullopt; // cannot use map anymore ! (it was moved)
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    auto source = trivis::geom::MakePoint(pov.x, pov.y);
    auto target = trivis::geom::MakePoint(pov.tx, pov.ty);
    auto direction = (target - source).CopyNormalized();
    auto vis_radius = pov.vis_radius > 0.0 ? std::make_optional(pov.vis_radius) : std::nullopt;

    LOGF_INF("Locating the source point " << source << ".");
    clock.Restart();
    auto source_pl = vis.LocatePoint(source);
    bool is_source_inside = source_pl.has_value();
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    bool visible = false;
    std::optional<trivis::geom::FPoint> ray_intersection;
    if (!is_source_inside) {
        LOGF_WRN("Source point " << source << " is outside of the map.");
    } else {
        if (pov.shoot_ray) {
            LOGF_INF("Shooting a ray from " << source << " in the direction of " << target << ".");
            auto ray_shooting_result = vis.ShootRay(source, source_pl.value(), direction);
            ray_intersection = ray_shooting_result.p;
            if (vis_radius && ray_intersection->DistanceTo(source) > vis_radius.value()) {
                visible = false;
            } else {
                visible = true;
            }
            LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
        } else {
            LOGF_INF("Determining visibility between " << source << " and " << target << ".");
            visible = vis.IsVisible(source, source_pl.value(), target, vis_radius);
            LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
        }
    }

    LOGF_INF("Preparing to draw the result.");
    clock.Restart();
    std::string pdf_file_str = pov.output_dir + "/";
    pdf_file_str += "eg_vis";
    if (pov.shoot_ray) {
        pdf_file_str += "_ray_shoot";
    } else {
        pdf_file_str += "_2point";
    }
    pdf_file_str += "_" + pov.map_name;
    pdf_file_str += "_" + source.ToString("", "-", "");
    pdf_file_str += "_" + target.ToString("", "-", "");
    if (vis_radius) {
        pdf_file_str += "_d-" + std::to_string(vis_radius.value());
    }
    pdf_file_str += ".pdf";
    fs::create_directories(pov.output_dir);
    namespace dr = trivis_plus::drawing;
    auto drawer = dr::MakeMapDrawer(vis.map());
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Drawing the result to " << pdf_file_str << ".");
    clock.Restart();
    drawer.OpenPDF(pdf_file_str);
    drawer.DrawMap();
    if (is_source_inside) {
        auto color_line = visible ? dr::kColorLimeGreen : dr::kColorRed;
        if (ray_intersection) {
            drawer.DrawLine(source, ray_intersection.value(), 0.2, color_line, 0.5);
        } else {
            drawer.DrawLine(source, target, 0.2, color_line, 0.5);
        }
    }
    if (vis_radius) {
        drawer.DrawArc(source, vis_radius.value(), 0.0, 2 * M_PI, 0.2, dr::kColorBlueViolet, 0.2);
    }
    drawer.DrawLine(source, source + direction, 0.2, dr::kColorBlueViolet);
    drawer.DrawPoint(source, 0.3, dr::kColorBlueViolet);
    double target_opacity = pov.shoot_ray ? 0.2 : 1.0;
    drawer.DrawPoint(target, 0.3, dr::kColorBlueViolet, target_opacity);
    if (ray_intersection) {
        drawer.DrawPoint(ray_intersection.value(), 0.3, dr::kColorRed);
    }
    drawer.Close();
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Finished.");

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