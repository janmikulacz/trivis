# TřiVis

***Versatile, Reliable, and High-Performance Tool for Computing Visibility in Polygonal Environments***

## About

TřiVis, named after the Czech word *'tři'* (meaning *'three'*) and the term *'visibility'*, is a C++ library for computing visibility-related structures in 2D polygonal environments.
It is based on the *triangular expansion algorithm* ([Bungiu et al., 2014](https://arxiv.org/abs/1403.3905); [Xu and Güting, 2015](https://doi.org/10.1007/s10707-014-0213-7)), which uses preprocessing to convert the input polygonal environment into a triangular mesh, and then traverses the mesh to compute visibility regions, visibility graphs, ray-shooting queries, and other visibility-related structures.

## Definitions

- **Polygonal Environment:** A 2D environment consisting of a single outer boundary and zero or more inner boundaries (holes), where each boundary is a simple polygon.
- **Triangular Mesh:** A mesh of triangles including adjacency information that covers the polygonal environment, where the intersection of any two triangles is either empty, a common edge, or a common vertex.
- **Point Location Query:** Determines the triangle of the triangular mesh that contains a given point in the polygonal environment or reports that the point is outside the environment.
- **Visibility:** Two points in the polygonal environment are visible to each other if the line segment connecting them lies entirely within the environment (touching the boundaries is allowed).
- **Two-Point Visibility Query:** Determines whether two given points in the polygonal environment are visible to each other.
- **Ray-Shooting Query:** Determines an intersection point of a ray starting at a given interior point in a given direction with the boundary of the polygonal environment.
- **Visibility Region Query:** Determines the region of the polygonal environment visible from a given interior point.
- **Visible Vertices Query:** Determines the subset of vertices of the polygonal environment visible from a given interior point.
- **Visible Points Query:** Determines the subset of given interior points that are visible from another given interior point.
- **Vertex-Vertex Visibility Graph:** A graph where vertices represent the vertices of the polygonal environment, and edges represent visibility between pairs of vertices.
- **Point-Point Visibility Graph:** A graph where vertices represent given interior points, and edges represent visibility between pairs of points.
- **Point-Vertex Visibility Graph:** A graph where vertices represent both given interior points and vertices of the polygonal environment, and edges represent visibility between pairs of points and vertices.
- **Limited Visibility Range:** A constraint on the visibility queries that restricts the visibility range to a given distance from the query point.

## Main Features

### Supported Queries

TřiVis supports all the queries defined above, each with a specialized, easy-to-use API.

### Worst-Case Query Complexity

The main feature of TřiVis is the *triangular expansion algorithm* ([Bungiu et al., 2014](https://arxiv.org/abs/1403.3905); [Xu and Güting, 2015](https://doi.org/10.1007/s10707-014-0213-7)), which provides $O(n)$ query time for two-point visibility and ray-shooting queries, $O(nh)$ for visibility region and visible vertices/points queries, $O(n^2h)$ for the vertex-vertex visibility graph, and $O(pnh)$ for point-point and point-vertex visibility graphs. 
Here, $n$ is the number of vertices in the polygonal environment, $h$ is the number of holes, and $p$ is the number of input interior points. 

It is important to note that the query time complexity is worst-case, and the actual query time may be significantly lower in practice. This is especially true for TřiVis, as the *triangular expansion algorithm* is known to exhibit output-sensitive behavior to some extent, traversing only the triangles that are partially or fully visible from the query point ([Bungiu et al., 2014](https://arxiv.org/abs/1403.3905)).

### Efficient Limited-Range Queries

TřiVis supports limited-range queries, employing an early termination mechanism to avoid traversing unnecessary triangles of the triangular mesh.
This can significantly reduce the query time for queries with a highly restricted visibility range relative to the average size of the open space in the environment.

### Fast Point Location

TřiVis employs a bucketing technique to accelerate point location queries, which are essential for visibility-related queries.
Using the bucketing technique from ([Edahiro et al., 1984](https://doi.org/10.1145/357337.357338)), TřiVis achieves $O(1)$ expected query time for point location queries.

### Robustness and Reliability

TřiVis relies on floating-point arithmetic and incorporates a unique combination of adaptive robust geometry predicates ([Shewchuk, 1997](https://doi.org/10.1007/PL00009321)) and $\epsilon$-geometry ([Fortune and Van Wyk, 1993](https://doi.org/10.1145/160985.161015)).
Thanks to this combination, TřiVis is characterized by reliable and predictable behavior, free from crashing or infinite looping, providing consistent outputs that align with user expectations, as supported by the results of our evaluation on exceptionally complex query instances.
For more information, refer to the [Performance Evaluation](#performance-evaluation) section.

### Performance on Large Instances

TřiVis has been tested on 34 complex polygonal environments with up to 8,320 vertices and 679 holes by computing visibility regions for 1,000 interior points in each environment.
With an average query time of 9±6μs, TřiVis outperforms other tested implementations by at least an order of magnitude, while keeping the preprocessing time below 20ms for the tested polygonal environments.
For more information, refer to the [Performance Evaluation](#performance-evaluation) section.

## Building and Linking the Library

TřiVis is built with [CMake](https://cmake.org/).
To make the library available to your C++ project, you can copy the source code in the [trivis](trivis) directory and include it in your CMake project.
```CMake
# In your CMakeLists.txt:
add_subdirectory(path_to_trivis) # Add the library to your project.
target_link_libraries(your_target PUBLIC Trivis) # Link the library to your target.
```

The build has been tested with the [GCC](https://gcc.gnu.org/) 12.3.0 compiler and [CMake](https://cmake.org/) 3.28.1 on the [Ubuntu](https://ubuntu.com/) 20.04.6 LTS operating system.
For your convenience, the [conda/trivis.yml](conda/trivis.yml) file describes [Conda](https://docs.conda.io/en/latest/) environment with the exact compiler and CMake versions used by the authors.

To create and test the environment, run the following in the root directory of this repository:
```bash
# Assuming you have Conda installed:
conda env create -f conda/trivis.yml # Create the environment.
conda activate trivis # Activate the environment.
gcc --version # Check the compiler version (12.3.0).
cmake --version # Check the CMake version (3.28.1).
mkdir build # Create the build directory.
cd build # Go to the build directory.
cmake ../trivis # Configure the build (ignore the warnings).
make # Build the library.
```

Apart from TřiVis's source code, this repository also contains TřiVis+, an extension adding support for data loading and visualization, and an example project depending on TřiVis+, demonstrating the library's usage.

For more information, refer to the respective [TřiVis+](#třivis) and [Examples](#examples) sections. 

## Basic Usage

### Initialization

To use the library in your C++ project, include the main header file:
```C++
#include "trivis/trivis.h"
```
TřiVis's main API is the `trivis::Trivis` class, which provides methods for visibility-related queries.
It is initialized with a polygonal environment represented by a `trivis::geom::PolyMap` object.
```C++
trivis::Trivis vis;
{   // Scope for the environment.
    trivis::geom::PolyMap environment;
    // TODO: Fill the environment.
    vis = trivis::Trivis{std::move(environment)};
}   // The environment was moved to vis.
```
The triangular mesh of the environment has just been constructed.
To inspect it, you can access the mesh directly:
```C++
const trivis::mesh::TriMesh &mesh = vis.mesh();
```

### Computing Visibility Queries

To perform visibility-related queries, you need to provide the query point, the point location result, and the visibility range (if limited).
```C++
trivis::geom::FPoint query;
// TODO: Fill the query.
std::optional<trivis::Trivis::PointLocationResult> plr = vis.LocatePoint(query);
if (plr.has_value()) {
    // The query point is inside the environment.
} else {
    // The query point is outside the environment.
}
std::optional<double> range;
// TODO: Set the visibility range or leave it unlimited.
```
**Two-point visibility query:**
```C++
trivis::geom::FPoint target;
// TODO: Fill the target.
if (plr.has_value()) {
    bool is_visible = vis.IsVisible(query, plr.value(), target, range);
}
```
**Ray-shooting query:**
```C++
trivis::geom::FPoint direction;
// TODO: Fill the direction.
if (plr.has_value()) {
    std::optional<trivis::Trivis::RayShootingResult> rsr = vis.ShootRay(query, plr.value(), direction, range);
    if (rsr.has_value()) {
        trivis::geom::FPoint intersection = rsr->p;
    } else {
        // The intersection point is outside the visibility range.
    }
}
```
**Visibility region query:**
```C++
if (plr.has_value()) {
    trivis::AbstractVisibilityRegion avr = vis.VisibilityRegion(query, plr.value(), range);
    trivis::RadialVisibilityRegion rvr = vis.ToRadialVisibilityRegion(avr);
    if (range.has_value()) {
        rvr.IntersectWithCircleCenteredAtSeed();
    }
    // Optional post-processing:
    if (range.has_value()) {
        vis_region.SampleArcEdges(M_PI / 180.0);
    }
    trivis::geom::FPolygon polygon_approx = vis_region->ToPolygon();
}
```
**Visible vertices query:**
```C++
if (plr.has_value()) {
    std::vector<int> visible_vertices = vis.VisibleVertices(query, plr.value(), range);
}
```
**Visible points query:**
```C++
std::vector<trivis::geom::FPoint> input_points;
// TODO: Fill the input points.
std::vector<std::optional<trivis::Trivis::PointLocationResult>> input_plrs;
input_plrs.reserve(input_points.size());
for (const auto &p : input_points) {
    input_plrs.push_back(vis.LocatePoint(p));
}
if (plr.has_value()) {
    std::vector<int> visible_points = vis.VisiblePoints(query, plr.value(), input_points, input_plrs, range);
}
```
**Visibility graphs:**
```C++
std::vector<std::vector<int>> vertex_vertex_graph = vis.VertexVertexVisibilityGraph(range);
std::vector<std::vector<int>> point_point_graph = vis.PointPointVisibilityGraph(input_points, input_plrs, range);
std::vector<std::vector<int>> point_vertex_graph = vis.PointVertexVisibilityGraph(input_points, input_plrs, range);
```

## Structure and Dependencies

## TřiVis+

## Examples

```bash
# Play with the 2-point/ray-shooting example:
./build-Release/examples/vis_2point 
./build-Release/examples/vis_2point --vis-radius 10
./build-Release/examples/vis_2point --shoot-ray
./build-Release/examples/vis_2point --help # to see all options
# Play with the visibility region example:
./build-Release/examples/vis_region
./build-Release/examples/vis_region --vis-radius 8
./build-Release/examples/vis_region --vis-radius 8 --sample-arc-edges 0.1
./build-Release/examples/vis_region --vis-radius 8 --sample-arc-edges 0.1 --to-polygon
./build-Release/examples/vis_region --help # to see all options
# Play with the visible vertices example:
./build-Release/examples/vis_vertices
./build-Release/examples/vis_vertices --reflex-only
./build-Release/examples/vis_vertices --vis-radius 6
./build-Release/examples/vis_vertices --help # to see all options
# Play with the visibility graph example:
./build-Release/examples/vis_graph
./build-Release/examples/vis_graph --reflex-only
./build-Release/examples/vis_graph --reflex-only --vis-radius 6
./build-Release/examples/vis_graph --point-vertex --reflex-only --vis-radius 6
./build-Release/examples/vis_graph --point-point
./build-Release/examples/vis_graph --help # to see all options
```

## Performance Evaluation

```bash
conda create env -f conda/trivis_with_boost_cairo_cgal.yml # unless you have it already
conda activate trivis # unless already activated
export TRIVIS_BUILD_TRIVIS_PLUS=1 # to build the dependency
export TRIVIS_BUILD_PERFORMANCE=1 # to build the performance tests
bash build.bash # to run the build
# creates the build-Release directory with the executables
# todo: copy the map instances to data/maps
cd performance # go to the performance directory
bash gen_points_all.bash # to generate the test points
# generates the test points in the ../data/points directory (also saves the CDT mesh for each map in the ../data/meshes directory)
bash test_all.bash # to construct the visibility regions
# generates the raw outputs in the outputs directory 
bash eval_all.bash # to compare the results with the reference
# generates CSV files with the evaluation results in the results/csv directory
rm outputs/* # you can now remove the raw outputs
cd results # go to the results directory
conda create env -f conda/datatable.yml # unless you have it already
conda activate datatable
python process1.py # to process the results
# generates the results.csv file, you can now inspect it
cd ../.. # go back to the root directory
conda deactivate # to deactivate the datatable environment
conda deactivate # to deactivate the trivis environment
```

## Documentation

## Version History

## Planned Development

## License

## Troubleshooting
