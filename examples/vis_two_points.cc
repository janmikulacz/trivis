/**
 * File:   draw_visibility_region.cc
 *
 * Date:   04.11.2022
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include <iomanip>

// === BOOST INCLUDES ===

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// === TRIVIS INCLUDES ===

#include "trivis/trivis.h"
#include "trivis/geom/generic_geom_types.h"
#include "trivis/geom/generic_geom_utils.h"
#include "trivis/utils/random_points.h"
#include "trivis/utils/simple_clock.h"

// === THIS PROJECT INCLUDES ===

#include "trivis_plus/data_loading/load_map.h"
#include "trivis_plus/drawing/fancy_drawing.h"
#include "trivis_plus/utils/log.h"

#ifndef DEFAULT_MAP_DIR
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_OUT_DIR
#define DEFAULT_OUT_DIR "."
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
    std::string map_name = "jari-huge";
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string output_dir = DEFAULT_OUT_DIR;
    double map_scale = 0.01;
    double bucket_size = 1.0;
    double x1 = 8.5;
    double y1 = 22.0;
    double x2 = 17.0;
    double y2 = 2.0;
    double vis_radius = -1.0;
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
         po::value(&pov.map_scale),
         "Map scale (optional).")
        ("bucket-size",
         po::value(&pov.bucket_size)->default_value(pov.bucket_size),
         "Size of the buckets used for locating a point in the triangular mesh.")
        ("x1",
         po::value(&pov.x1)->default_value(pov.x1),
         "X coordinate of the first query point.")
        ("y1",
         po::value(&pov.y1)->default_value(pov.y1),
         "Y coordinate of the first query point.")
        ("x2",
         po::value(&pov.x2)->default_value(pov.x2),
         "X coordinate of the second query point.")
        ("y2",
         po::value(&pov.y2)->default_value(pov.y2),
         "Y coordinate of the second query point.")
        ("vis-radius",
         po::value(&pov.vis_radius)->default_value(pov.vis_radius),
         "Limited visibility radius (-1 ~ infinite).");
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

    LOGF_INF("HELLO WORLD!");

    fs::create_directories(pov.output_dir);

    trivis::utils::SimpleClock clock;

    LOGF_INF(">> Initializing Trivis.");
    clock.Restart();
    // Create and initialize the Trivis object.
    trivis::Trivis vis;
    {   // Load map from file and move it to Trivis (without copying).
        trivis::geom::PolyMap map;
        std::string load_msg = trivis_plus::data_loading::LoadPolyMapSafely(pov.map_dir + "/" + pov.map_name + pov.map_extension, map, pov.map_scale);
        if (load_msg != "ok") {
            LOGF_FTL("Error while loading map. " << load_msg);
            return EXIT_FAILURE;
        }
        // === Bellow is some optional map preprocessing. ===
        // Warning: order of the following matters!
        map.ShiftToOrigin(); // Subtracts min X and Y coordinates from all points.
        map.RemoveDuplicatePoints(); // Removes all consecutive identical points.
        map.SimplifyWeaklySimplePolygons(); // Splits all weakly simple polygons to multiple (touching) strongly simple polygons.
        map.RemoveCollinearPoints(); // Removes all consecutive collinear points.
        // ==================================================
        vis.SetMap(std::move(map)); // Set the map to Trivis instance.
        // cannot use map anymore ! (it was moved)
    }
    vis.ConstructMeshCDT();
    vis.FillPointLocationBuckets(pov.bucket_size);
    //vis.OptimizePointLocationBucketTriangles(); // Optional to slightly improve bucketing speed.
    LOGF_INF("<< DONE. It took " << clock.TimeInSeconds() << " seconds.");

    const auto &lim = vis.limits();
    LOGF_INF("Map limits [ MIN: " << lim.x_min << ", " << lim.y_min << " | MAX: " << lim.x_max << ", " << lim.y_max << " ].");

    // Determine visibility.
    clock.Restart();
    auto vis_radius_opt = pov.vis_radius > 0.0 ? std::make_optional(pov.vis_radius) : std::nullopt;
    auto q = trivis::geom::MakePoint(pov.x1, pov.y1);
    auto p = trivis::geom::MakePoint(pov.x2, pov.y2);
    LOGF_INF(">> Computing visibility between " << q << " and " << p << ".");
    auto q_location_opt = vis.LocatePoint(q);
    if (!q_location_opt) {
        LOGF_FTL("Point " << q << " is outside of the map.");
        return EXIT_FAILURE;
    }
    std::optional<bool> visible_opt = vis.IsVisible(q, *q_location_opt, p, vis_radius_opt);
    if (!visible_opt) {
        LOGF_FTL("Error while computing visibility.");
        return EXIT_FAILURE;
    }
    bool visible = visible_opt.value_or(false);
    LOGF_INF("<< DONE. It took " << clock.TimeInSeconds() << " seconds.");

    // Draw the result.
    std::string pdf_file_str = pov.output_dir + "/";
    pdf_file_str += "eg_vis_2points";
    pdf_file_str += "_" + pov.map_name;
    pdf_file_str += "_" + q.ToString("", "-", "");
    pdf_file_str += "_" + p.ToString("", "-", "");
    if (vis_radius_opt) {
        pdf_file_str += "_d-" + std::to_string(*vis_radius_opt);
    }
    pdf_file_str += ".pdf";
    LOGF_INF("Drawing the result to " << pdf_file_str << ".");
    auto drawer = dr::MakeMapDrawer(vis.map());
    drawer.OpenPDF(pdf_file_str);
    dr::FancyDrawMap(drawer, vis);
    dr::RGB color_segment = visible ? dr::kColorLimeGreen : dr::kColorRed;
    drawer.DrawLine(q, p, 0.1, color_segment);
    drawer.DrawPoint(q, 0.1, dr::kColorBlueViolet);
    if (vis_radius_opt) {
        drawer.DrawArc(q, *vis_radius_opt, 0.0, 2 * M_PI, 0.05, dr::kColorBlueViolet);
    }
    drawer.DrawPoint(p, 0.1, dr::kColorNavy);
    drawer.Close();

    LOGF_INF("BYE BYE WORLD!");

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