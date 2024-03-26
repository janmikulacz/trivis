/**
 * File:   sandbox.cc
 *
 * Date:   26.03.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */


#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "trivis/trivis.h"

#include "trivis_plus/trivis_plus.h"

#ifndef DEFAULT_MAP_DIR // May be defined in CMakeLists.txt.
#define DEFAULT_MAP_DIR "."
#endif

#ifndef DEFAULT_OUT_DIR // May be defined in CMakeLists.txt.
#define DEFAULT_OUT_DIR "."
#endif

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct ProgramOptionVariables {
    std::string verbosity = "info";
    std::string map_name = "potholes";
    std::string map_extension = ".txt";
    std::string map_dir = DEFAULT_MAP_DIR;
    std::string map_full_path;
    std::string out_dir = DEFAULT_OUT_DIR;
    std::string out_pdf;
    double map_scale = -1.0;
    double vis_radius = -1.0;
    // TODO: add more variables
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
         po::value(&pov.out_dir)->default_value(pov.out_dir),
         "Output directory.")
        ("out-pdf",
         po::value(&pov.out_pdf)->default_value(pov.out_pdf),
         "Output pdf file.")
        ("map-scale",
         po::value(&pov.map_scale)->default_value(pov.map_scale),
         "Map coordinates are scaled by this factor when loading the map (< 0.0 ~ no scaling or scale is loaded with the map).")
        ("vis-radius",
         po::value(&pov.vis_radius)->default_value(pov.vis_radius),
         "Limited visibility radius (-1 ~ infinite).")
        // TODO: add more options
        ;
}

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

char ParseProgramOptions(
    int argc,
    const char *const *argv,
    ProgramOptionVariables &pov
) {
    po::variables_map vm;
    po::options_description command_line_options;
    po::options_description options_description("Program options");
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

int MainBody(const ProgramOptionVariables &pov) {

    LOGF_INF("Running the visibility region example.");

    trivis::utils::SimpleClock clock;
    std::stringstream info;

    LOGF_INF("Loading map from " << pov.map_full_path << ".");
    clock.Restart();
    std::optional<double> scale = pov.map_scale > 0.0 ? std::make_optional(pov.map_scale) : std::nullopt;
    auto map = trivis_plus::data_loading::LoadPolyMap(pov.map_full_path, scale, &info);
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

    const auto &lim = vis.limits();
    LOGF_INF("Map limits: [ MIN: " << lim.x_min << ", " << lim.y_min << " | MAX: " << lim.x_max << ", " << lim.y_max << " ].");

    // TODO: do whatever you want

    LOGF_INF("Preparing to draw the result.");
    clock.Restart();
    std::string pdf_file_str = pov.out_dir + "/";
    if (pov.out_pdf.empty()) {
        pdf_file_str += "sandbox.pdf";
    } else {
        pdf_file_str += pov.out_pdf;
    }
    fs::create_directories(pov.out_dir);
    namespace dr = trivis_plus::drawing;
    auto drawer = dr::MakeMapDrawer(vis.map());
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    LOGF_INF("Drawing the result to " << pdf_file_str << ".");
    clock.Restart();
    drawer.OpenPDF(pdf_file_str);
    dr::FancyDrawMap(drawer, vis);
    // TODO: draw whatever you want
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