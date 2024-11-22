/**
 * File:   pathfinder.cc
 *
 * Date:   21.11.2024
 * Author: Jan Mikula
 * E-mail: jan.mikula@cvut.cz
 *
 */

#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "trivis/trivis.h"

#include "trivis_plus/trivis_plus.h"

#include "trivis_pathfinder/trivis_pathfinder.h"

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
    double x = 4.0;
    double y = 0.5;
    double tx = 19.0;
    double ty = 19.0;
    bool all_pairs_reflex = false;
    bool all_pairs_random_cities = false;
    bool all_pairs_random_points = false;
    bool precompute_reflex_paths = false;
    bool precompute_cities_graph = false;
    bool precompute_cities_paths = false;
    int n_rand_points = 20;
    int rand_seed = 0;
    bool lengths_only = false;
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
         "Map coordinates are scaled by this factor when loading the map (<= 0.0 ~ no scaling or scale is loaded with the map).")
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
        ("all-pairs-reflex",
         po::bool_switch(&pov.all_pairs_reflex)->default_value(pov.all_pairs_reflex),
         "Compute all pairs shortest paths between reflex vertices.")
        ("all-pairs-random-cities",
         po::bool_switch(&pov.all_pairs_random_cities)->default_value(pov.all_pairs_random_cities),
         "Compute all pairs shortest paths between random cities.")
        ("all-pairs-random-points",
         po::bool_switch(&pov.all_pairs_random_points)->default_value(pov.all_pairs_random_points),
         "Compute all pairs shortest paths between random points.")
        ("precompute-reflex-paths",
         po::bool_switch(&pov.precompute_reflex_paths)->default_value(pov.precompute_reflex_paths),
         "Precompute reflex shortest paths.")
        ("precompute-cities-graph",
         po::bool_switch(&pov.precompute_cities_graph)->default_value(pov.precompute_cities_graph),
         "Precompute cities visibility graph.")
        ("precompute-cities-paths",
         po::bool_switch(&pov.precompute_cities_paths)->default_value(pov.precompute_cities_paths),
         "Precompute cities shortest paths.")
        ("n-rand-points",
         po::value(&pov.n_rand_points)->default_value(pov.n_rand_points),
         "Number of random points.")
        ("rand-seed",
         po::value(&pov.rand_seed)->default_value(pov.rand_seed),
         "Random seed.")
        ("lengths-only",
         po::bool_switch(&pov.lengths_only)->default_value(pov.lengths_only),
         "Compute only lengths of the shortest paths.");
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

struct ShortestPathWithInfo {
    double length = -1.0;
    std::optional<int> id_source;
    std::optional<int> id_target;
    trivis::geom::FPoint source;
    trivis::geom::FPoint target;
    std::vector<int> id_path_no_endpoints;
    trivis::geom::FPoints points_path;
};

int MainBody(const ProgramOptionVariables &pov) {

    LOGF_INF("Running the visibility region example.");

    if (pov.all_pairs_reflex && pov.all_pairs_random_cities) {
        LOGF_FTL("Options --all-pairs-reflex and --all-pairs-random-cities cannot be used together.");
        return EXIT_FAILURE;
    }
    if (pov.all_pairs_reflex && pov.all_pairs_random_points) {
        LOGF_FTL("Options --all-pairs-reflex and --all-pairs-random-points cannot be used together.");
        return EXIT_FAILURE;
    }
    if (pov.all_pairs_random_cities && pov.all_pairs_random_points) {
        LOGF_FTL("Options --all-pairs-random-cities and --all-pairs-random-points cannot be used together.");
        return EXIT_FAILURE;
    }

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

    std::optional<trivis::geom::FPoint> source;
    std::optional<trivis::geom::FPoint> target;
    std::optional<trivis::geom::FPoints> query_points;
    if (pov.all_pairs_reflex) {
        // do nothing
    } else if (pov.all_pairs_random_cities || pov.all_pairs_random_points) {
        LOGF_INF("Generating " << pov.n_rand_points << " random query points.");
        clock.Restart();
        std::vector<double> accum_areas;
        auto rng = std::mt19937(pov.rand_seed);
        query_points = trivis::geom::FPoints{};
        query_points->reserve(pov.n_rand_points);
        for (int i = 0; i < pov.n_rand_points; ++i) {
            query_points->push_back(trivis::utils::UniformRandomPointInRandomTriangle(vis.triangles(), accum_areas, rng));
        }
        LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
    } else {
        source = std::make_optional(trivis::geom::FPoint(pov.x, pov.y));
        target = std::make_optional(trivis::geom::FPoint(pov.tx, pov.ty));
    }

    trivis::utils::SimpleClock clock_total;

    LOGF_INF("Constructing reflex visibility graph.");
    clock.Restart();
    trivis_pathfinder::TrivisPathfinder pathfinder;
    auto status = pathfinder.ConstructReflexVisibilityGraph(vis);
    if (status != trivis_pathfinder::utils::Status::kOk) {
        LOGF_FTL("Error while constructing reflex visibility graph.");
        return EXIT_FAILURE;
    }
    LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    if (pov.precompute_reflex_paths) {
        LOGF_INF("Precomputing reflex shortest paths.");
        clock.Restart();
        status = pathfinder.PrecomputeReflexShortestPaths();
        if (status != trivis_pathfinder::utils::Status::kOk) {
            LOGF_FTL("Error while precomputing reflex shortest paths.");
            return EXIT_FAILURE;
        }
        LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
    }

    std::vector<std::vector<double>> lengths;
    std::vector<ShortestPathWithInfo> shortest_paths;
    if (pov.all_pairs_reflex) {
        LOGF_INF("Computing all pairs shortest paths between reflex vertices.");
        clock.Restart();
        if (pov.lengths_only) {
            lengths = std::vector<std::vector<double>>(pathfinder.n_reflex(), std::vector<double>(pathfinder.n_reflex(), 0.0));
            for (int i = 0; i < pathfinder.n_reflex(); ++i) {
                for (int j = i + 1; j < pathfinder.n_reflex(); ++j) {
                    trivis_pathfinder::utils::StatusWithResult<double> result;
                    result = pathfinder.ShortestPathReflex(vis, i, j);
                    if (result.status != trivis_pathfinder::utils::Status::kOk) {
                        LOGF_FTL("Error while computing the shortest path from reflex vertex " << i << " to reflex vertex " << j << ".");
                        return EXIT_FAILURE;
                    }
                    lengths[i][j] = result.result;
                    lengths[j][i] = result.result;
                }
            }
        } else {
            shortest_paths.reserve(pathfinder.n_reflex() * (pathfinder.n_reflex() - 1) / 2);
            for (int i = 0; i < pathfinder.n_reflex(); ++i) {
                for (int j = i + 1; j < pathfinder.n_reflex(); ++j) {
                    ShortestPathWithInfo sp;
                    sp.id_source = i;
                    sp.id_target = j;
                    sp.source = pathfinder.reflex_points()[i];
                    sp.target = pathfinder.reflex_points()[j];
                    auto result = pathfinder.ShortestPathReflex(vis, i, j, &sp.points_path, &sp.id_path_no_endpoints);
                    if (result.status != trivis_pathfinder::utils::Status::kOk) {
                        LOGF_FTL("Error while computing the shortest path from reflex vertex " << i << " to reflex vertex " << j << ".");
                        return EXIT_FAILURE;
                    }
                    sp.length = result.result;
                    shortest_paths.push_back(std::move(sp));
                }
            }
        }
        LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    } else if (pov.all_pairs_random_cities || pov.all_pairs_random_points) {

        if (pov.all_pairs_random_cities || pov.precompute_cities_graph) {
            LOGF_INF("Constructing cities visibility graph.");
            clock.Restart();
            status = pathfinder.ConstructCitiesVisibilityGraph(vis, query_points.value());
            if (status != trivis_pathfinder::utils::Status::kOk) {
                LOGF_FTL("Error while constructing cities visibility graph.");
                return EXIT_FAILURE;
            }
            LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
        }

        if (pov.precompute_cities_paths) {
            LOGF_INF("Precomputing cities shortest paths.");
            clock.Restart();
            std::vector<std::vector<int>> id_paths_no_endpoints;
            std::vector<double> paths_lengths;
            status = pathfinder.PrecomputeCitiesShortestPaths();
            if (status != trivis_pathfinder::utils::Status::kOk) {
                LOGF_FTL("Error while precomputing cities shortest paths.");
                return EXIT_FAILURE;
            }
            LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
        }

        LOGF_INF("Computing all pairs shortest paths between random points.");
        clock.Restart();
        if (pov.lengths_only) {
            lengths = std::vector<std::vector<double>>(query_points->size(), std::vector<double>(query_points->size(), 0.0));
            for (int i = 0; i < query_points->size(); ++i) {
                const auto &point_i = query_points->operator[](i);
                for (int j = i + 1; j < query_points->size(); ++j) {
                    const auto &point_j = query_points->operator[](j);
                    trivis_pathfinder::utils::StatusWithResult<double> result;
                    if (pov.all_pairs_random_cities) {
                        result = pathfinder.ShortestPathCities(vis, i, j);
                    } else {
                        result = pathfinder.ShortestPathPoints(vis, point_i, point_j);
                    }
                    if (result.status != trivis_pathfinder::utils::Status::kOk) {
                        LOGF_FTL("Error while computing the shortest path from random point " << i << " to random point " << j << ".");
                        return EXIT_FAILURE;
                    }
                    lengths[i][j] = result.result;
                    lengths[j][i] = result.result;
                }
            }
        } else {
            shortest_paths.reserve(query_points->size() * (query_points->size() - 1) / 2);
            for (int i = 0; i < query_points->size(); ++i) {
                const auto &point_i = query_points->operator[](i);
                for (int j = i + 1; j < query_points->size(); ++j) {
                    const auto &point_j = query_points->operator[](j);
                    ShortestPathWithInfo sp;
                    sp.id_source = i;
                    sp.id_target = j;
                    sp.source = point_i;
                    sp.target = point_j;
                    trivis_pathfinder::utils::StatusWithResult<double> result;
                    if (pov.all_pairs_random_cities) {
                        result = pathfinder.ShortestPathCities(vis, i, j, &sp.points_path, &sp.id_path_no_endpoints);
                    } else {
                        result = pathfinder.ShortestPathPoints(vis, point_i, point_j, &sp.points_path, &sp.id_path_no_endpoints);
                    }
                    if (result.status != trivis_pathfinder::utils::Status::kOk) {
                        LOGF_FTL("Error while computing the shortest path from random point " << i << " to random point " << j << ".");
                        return EXIT_FAILURE;
                    }
                    sp.length = result.result;
                    shortest_paths.push_back(std::move(sp));
                }
            }
        }
        LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    } else {

        LOGF_INF("Computing the shortest path.");
        clock.Restart();
        if (pov.lengths_only) {
            lengths = std::vector<std::vector<double>>(2, std::vector<double>(2, 0.0));
            trivis_pathfinder::utils::StatusWithResult<double> result;
            result = pathfinder.ShortestPathPoints(vis, source.value(), target.value());
            if (result.status != trivis_pathfinder::utils::Status::kOk) {
                LOGF_FTL("Error while computing the shortest path.");
                return EXIT_FAILURE;
            }
            lengths[0][1] = result.result;
            lengths[1][0] = result.result;
        } else {
            shortest_paths.reserve(1);
            ShortestPathWithInfo sp;
            sp.source = source.value();
            sp.target = target.value();
            auto result = pathfinder.ShortestPathPoints(vis, source.value(), target.value(), &sp.points_path, &sp.id_path_no_endpoints);
            if (result.status != trivis_pathfinder::utils::Status::kOk) {
                LOGF_FTL("Error while computing the shortest path.");
                return EXIT_FAILURE;
            }
            sp.length = result.result;
            shortest_paths.push_back(sp);
        }
        LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");

    }

    LOGF_INF("Total time: " << clock_total.TimeInSeconds() << " seconds.");

    if (pov.lengths_only) {
        LOGF_INF("Nothing to draw because --lengths-only option is set.");
    } else {
        LOGF_INF("Preparing to draw the result.");
        clock.Restart();
        std::string pdf_file_str = pov.out_dir + "/";
        if (pov.out_pdf.empty()) {
            pdf_file_str += "ex_pathfinder";
            if (pov.all_pairs_reflex) {
                pdf_file_str += "_ap_reflex"; // ap ~ all_pairs
            } else if (pov.all_pairs_random_cities) {
                pdf_file_str += "_ap_rand_cities"; // ap ~ all_pairs
            } else if (pov.all_pairs_random_points) {
                pdf_file_str += "_ap_rand_points"; // ap ~ all_pairs
            }
            pdf_file_str += "_" + pov.map_name;
            if (source.has_value()) {
                pdf_file_str += "_" + source->ToString("", "-", "");
            }
            if (target.has_value()) {
                pdf_file_str += "_" + target->ToString("", "-", "");
            }
            if (query_points.has_value()) {
                pdf_file_str += "_n-" + std::to_string(pov.n_rand_points);
                pdf_file_str += "_s-" + std::to_string(pov.rand_seed);
            }
            if (pov.precompute_reflex_paths) {
                pdf_file_str += "_prp"; // prp ~ precompute_reflex_paths
            }
            if (query_points.has_value()) {
                if (pov.all_pairs_random_points) {
                    if (pov.precompute_cities_graph) {
                        pdf_file_str += "_pcg"; // pcg ~ precompute_cities_graph
                    }
                }
                if (pov.precompute_cities_paths) {
                    pdf_file_str += "_pcp"; // pcp ~ precompute_cities_paths
                }
            }
            pdf_file_str += ".pdf";
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
        drawer.DrawMap();
        if (shortest_paths.size() == 1) {
            drawer.DrawPath(shortest_paths.front().points_path, 0.2, dr::kColorLimeGreen, 0.5);
        } else {
            for (const auto &sp: shortest_paths) {
                drawer.DrawPath(sp.points_path, 0.05, dr::kColorLimeGreen, 0.25);
            }
        }
        drawer.Close();
        LOGF_INF("DONE. It took " << clock.TimeInSeconds() << " seconds.");
    }

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